(ns mgs-bias.core
  (:require
   [cljs.pprint :as pp]
   [goog.string :as gstring]
   [goog.string.format]
   [goog.dom :as gdom]
   [reagent.core :as r]
   [cljsjs.d3]))

(defn say-hi []
  "Hi from core.cljs!")

(def html-id-bias-table "mgs-bias-protocol-steps")
(def html-id-num-otus "mgs-bias-num-otus")
(def html-id-otu-table-actual-counts "mgs-bias-otu-table-actual-counts")
(def html-id-otu-table-observed-counts "mgs-bias-otu-table-observed-counts")
(def html-id-composition-charts "mgs-bias-composition-charts")
(def html-id-composition-charts-svg "mgs-bias-composition-charts-svg")
(def html-id-update-charts-button "mgs-bias-update-charts")

(def config (atom {:num-otus 3}))

;;;; Data

(def otus [:otu-1 :otu-2 :otu-3])

(def protocol-steps [:extraction-bais
                     :pcr-bias
                     :sequencing-bias
                     :bioinformatics-bais])

;; row i is OTU i, col j is sample j
(def actual-otu-counts
  (r/atom [[1 15] ; OTU 1
           [1  1] ; OTU 2
           [1  4] ; OTU 3
           ]))

(def observed-otu-counts
  (r/atom nil))

;; row i is OTU i, col j is protocol step j
(def protocol-otu-bias-per-step
  (r/atom [[ 1 1 1    1]   ; OTU 1
           [ 4 6 0.5  1.5] ; OTU 2
           [15 2 0.25 0.8] ; OTU 3
           ]))

(def num-otus 3)
(def num-samples 2)

;; (def protocol-otu-bias-total
;;   (r/atom nil))

(defn total-bias
  "Returns a vec with the total bias for each OTU."
  [protocol-bias]
  (vec (map #(reduce * %) protocol-bias)))

(defn observed-counts
  "Returns vec of vecs to match the input data format."
  [actual-counts protocol-bias]
  (vec (map (fn [[otu-bias per-sample-counts]]
                (vec (map #(* otu-bias %) per-sample-counts)))
              (zipmap (total-bias protocol-bias)
                      actual-counts))))

;;;; Utils

(defn event-val [event]
  (-> event .-target .-value))

(defn counts-to-rel-abund [counts]
  (vec (map #(/ % (reduce + counts)) counts)))

(defn transpose
  "Transpose a matrix (2d vector)"
  [mat]
  (apply map vector mat))

(defn calc-rel-abund
  "Given a 2d count vec with rows OTU, cols sample (entry ij is otu i
  count in sample j), return a vec of relative abundances (entry ij is
  otu i relative abundance in sample j)."
  [counts]
  (defn rabund [vec]
    (let [total (apply + vec)]
      (map #(/ % total) vec)))

  (let [;; samples => 2d vec, rows are samples, cols are OTUs
        samples (transpose counts)]
    (transpose (map rabund samples))))

(def total-bar-height 350)
(defn calc-rect-coords
  "rel-abunds is a vector of relative abundances for all the OTUs for
  a single sample.  Some OTUs could have 0 relative abundance.  All
  vals should run from 0 to 1."
  [rel-abunds]
  (let [heights (map #(* total-bar-height %) rel-abunds)]
    (map-indexed (fn [i n]
                   (if (zero? i)
                     {:start 0 :height n}
                     ;; Add up all the heights before this one
                     {:start (apply + (take i heights))
                      :height n}))
                 heights)))

(def otu-colors ["#004488" "#DDAA33" "#BB5566"])

(defn rabund-bar [sample-idx rect-coords]
  (map-indexed (fn [i coords]
                 (let [x (* 100 (inc sample-idx))
                       y (:start coords)
                       width 50
                       height (:height coords)]
                   ^{:key (str x "-" y)}
                   [:rect {:x x :y y
                           :width width :height height
                           :fill (otu-colors i)}]))
               rect-coords))

(defn bar-charts
  "Both arguments are atoms."
  [otu-counts otu-bias-per-step]
  (let [rel-abund-by-sample (transpose (calc-rel-abund @otu-counts))
        all-rect-coords (map calc-rect-coords rel-abund-by-sample)]

    (print "@otu-counts")
    (pp/pprint @otu-counts)

    (print "(calc-rel-abund @otu-counts)")
    (pp/pprint (calc-rel-abund @otu-counts))

    (print "rel-abund-by-sample")
    (pp/pprint rel-abund-by-sample)

    (print "all-rect-coords")
    (pp/pprint all-rect-coords)
    ;; (print (map rabund-bar (range num-samples) all-rect-coords))
    [:svg {:width 500 :height 500}
     (map rabund-bar (range num-samples) all-rect-coords)]))



#_(defn bar-charts [otu-counts otu-bias-per-step]
  ;; Make sure any existing svg has been removed
  (-> js/d3
      (.select (str "#" html-id-composition-charts-svg))
      .remove)
  (let [svg (-> js/d3
                (.select (str "#" html-id-composition-charts))
                (.append "svg")
                (.attr "id" html-id-composition-charts-svg)
                (.attr "width" 500)
                (.attr "height" 500))]
    (-> svg
        (.append "circle")
        (.attr "r" 10)
        (.attr "cx" 100)
        (.attr "cy" 100))))

;;;; Components

(defn bias-table-header []
  [:thead
   [:tr
    [:th "OTU ID"]
    [:th "Extraction bias"]
    [:th "PCR bias"]
    [:th "Sequencing bias"]
    [:th "Bioinformatics bias"]
    [:th "Total bias"]]])

(defn check-NaN [event]
  (let [val (js/parseFloat (event-val event))]
    (if (js/isNaN val)
      0
      val)))

(defn bias-table-body
  "Assumes all OTUs have equal number of steps."
  [otu-bias-per-step]
  (let [num-otus (count @otu-bias-per-step)
        num-steps (count (first @otu-bias-per-step))]
    [:tbody
     (doall
      (for [otu-idx (range num-otus)]
        ^{:key (str "otu-" otu-idx)}
        [:tr
         [:td (str "otu-" (inc otu-idx))]
         (doall
          (for [step-idx (range num-steps)]
            ^{:key (str "step-" step-idx)}
            [:td
             [:input {:type "text"
                      :value (get-in @otu-bias-per-step
                                     [otu-idx step-idx])
                      :on-change (fn [event]
                                   (let [val (check-NaN event)]
                                     (swap! otu-bias-per-step
                                            assoc-in
                                            [otu-idx step-idx]
                                            val)
                                     (reset! observed-otu-counts
                                             (observed-counts @actual-otu-counts
                                                              @protocol-otu-bias-per-step))))}]]))
         [:td (reduce * (get @otu-bias-per-step otu-idx))]]))]))

(defn bias-table [data]
  [:table
   [bias-table-header]
   [bias-table-body data]])

(defn count-table-header []
  [:thead
   [:tr
    [:th "OTU ID"]
    [:th "Sample 1"]
    [:th "sample 2"]]])

;; parseFloat -> need to handle NaN
(defn actual-count-table-body
  [otu-counts]
  (let [num-otus (count @otu-counts)
        num-samples (count (first @otu-counts))]
    [:tbody
     (doall
      (for [otu-idx (range num-otus)]
        ^{:key (str "otu-" otu-idx)}
        [:tr
         [:td (str "otu-" (inc otu-idx))]
         (doall
          (for [sample-idx (range num-samples)]
            ^{:key (str "sample-" sample-idx)}
            [:td
             [:input {:type "text"
                      :value (get-in @otu-counts [otu-idx sample-idx])
                      :on-change (fn [event]
                                   (let [val (check-NaN event)]
                                     (swap! otu-counts
                                            assoc-in
                                            [otu-idx sample-idx]
                                            val)
                                     (reset! observed-otu-counts
                                             (observed-counts @actual-otu-counts
                                                              @protocol-otu-bias-per-step))))}]]))]))]))

(defn actual-count-table [data]
  [:table
   [count-table-header]
   [actual-count-table-body data]])

#_(defn observed-count-table-body
  "Assumes the data is all in the proper shape."
  [otu-counts otu-bias-per-step]
  (let [num-otus (count @otu-counts)
        num-samples (count (first @otu-counts))
        bias (total-bias @otu-bias-per-step)]
    [:tbody
     (doall
      (for [otu-idx (range num-otus)]
        ^{:key (str "otu-" otu-idx)}
        [:tr
         [:td (str "otu-" (inc otu-idx))]
         (doall
          (for [sample-idx (range num-samples)]
            (let [count (get-in @otu-counts [otu-idx sample-idx])]
              ^{:key (str "sample-" sample-idx)}
              [:td (* count (get bias otu-idx))])))]))]))

(defn observed-count-table-body
  "Assumes the data is all in the proper shape."
  [otu-counts]
  [:tbody
     (doall
      (for [otu-idx (range num-otus)]
        ^{:key (str "otu-" otu-idx)}
        [:tr
         [:td (str "otu-" (inc otu-idx))]
         (doall
          (for [sample-idx (range num-samples)]
            (let [count (get-in @otu-counts [otu-idx sample-idx])]
              ^{:key (str "sample-" sample-idx)}
              [:td count])))]))])

(defn observed-count-table [otu-counts otu-bias-per-step]
  [:table
   [count-table-header]
   [observed-count-table-body otu-counts otu-bias-per-step]])

#_(defn update-charts []
  [:input {:type "button"
           :value "Update charts!"
           :on-click (fn [e]
                       (bar-charts actual-otu-counts
                                        protocol-otu-bias-per-step))}])

;;;; Rendering

(reset! observed-otu-counts
        (observed-counts @actual-otu-counts
                         @protocol-otu-bias-per-step))

(r/render [bias-table protocol-otu-bias-per-step]
          (.getElementById js/document html-id-bias-table))

(r/render [actual-count-table actual-otu-counts]
          (.getElementById js/document html-id-otu-table-actual-counts))

(r/render [observed-count-table
           observed-otu-counts]
          (.getElementById js/document
                           html-id-otu-table-observed-counts))

(r/render [bar-charts actual-otu-counts observed-otu-counts]
          (.getElementById js/document html-id-composition-charts))

;; (r/render [update-charts]
;;           (.getElementById js/document html-id-update-charts-button))
