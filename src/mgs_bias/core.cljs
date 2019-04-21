(ns mgs-bias.core
  (:require
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

(def config (atom {:num-otus 3}))

;;;; Data

(def otus [:otu-1 :otu-2 :otu-3])

(def protocol-steps [:extraction-bais
                     :pcr-bias
                     :sequencing-bias
                     :bioinformatics-bais])

;; row i is OTU i, col j is sample j
(def actual-otu-counts
  (r/atom [[1 15]
           [1 1]
           [1 4]]))

(def observed-otu-counts
  (r/atom nil))

;; row i is OTU i, col j is protocol step j
(def protocol-otu-bias-per-step
  (r/atom [[1  1 1    1]
           [4  6 0.5  1.5]
           [15 2 0.25 0.8]]))

(def protocol-otu-bias-total
  (r/atom nil))

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
                      :value (get-in @otu-bias-per-step [otu-idx step-idx])
                      :on-change #(swap! otu-bias-per-step
                                         assoc-in
                                         [otu-idx step-idx]
                                         (event-val %))}]]))
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

(defn count-table-body
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
                      :on-change #(swap! otu-counts
                                         assoc-in
                                         [otu-idx sample-idx]
                                         (event-val %))}]]))]))]))

(defn count-table [data]
  [:table
   [count-table-header]
   [count-table-body data]])

(defn actual-count-table-body
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

(defn actual-count-table [otu-counts otu-bias-per-step]
  [:table
   [count-table-header]
   [actual-count-table-body otu-counts otu-bias-per-step]])

;;;; Rendering

(r/render [bias-table protocol-otu-bias-per-step]
          (.getElementById js/document html-id-bias-table))

(r/render [count-table actual-otu-counts]
          (.getElementById js/document html-id-otu-table-actual-counts))

(r/render [actual-count-table actual-otu-counts protocol-otu-bias-per-step]
          (.getElementById js/document html-id-otu-table-observed-counts))
