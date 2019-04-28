(ns mgs-bias.core
  (:require
   [goog.string :as gstring]
   [goog.string.format]
   [reagent.core :as r]))

(defn say-hi []
  "Hi from core.cljs!")

(def html-id-app "mgs-bias-app")
(def html-id-bias-table "mgs-bias-protocol-steps")
(def html-id-num-otus "mgs-bias-num-otus")
(def html-id-otu-table-actual-counts "mgs-bias-otu-table-actual-counts")
(def html-id-otu-table-observed-counts "mgs-bias-otu-table-observed-counts")
(def html-id-composition-charts "mgs-bias-composition-charts")
(def html-id-composition-charts-svg "mgs-bias-composition-charts-svg")
(def html-id-beta-table "mgs-bias-beta-table")

(def config (atom {:num-otus 3}))

;;;; Statistics

(defn rabund [v]
  (let [total (apply + v)]
    (vec  (map #(/ % total) v))))

(defn aitchison-dist
  "Input is OTU counts by sample (two 1D vecs).  Should be the same length."
  [xs ys]
  (assert (= (count xs) (count ys)))
  ;; Shadow original vars with their relative abundance.  If inputs
  ;; are already rel abundance this won't change the answer, but will
  ;; waste some cycles.
  (let [xs (rabund xs)
        ys (rabund ys)]
    (let [size (count xs)]
      (* (/ 1
            (* 2 size))
         (apply +
                (for [i (range size)
                      j (range size)]
                  (Math/pow (- (Math/log (/ (xs i)
                                            (xs j)))
                               (Math/log (/ (ys i)
                                            (ys j))))
                            2)))))))

(defn euclidean-dist
  "Same input as adist"
  [xs ys]
  (assert (= (count xs) (count ys)))
  (Math/sqrt
   (apply +
          (map-indexed (fn [i x]
                         (let [y (ys i)]
                           (Math/pow (- x y) 2)))
                       xs))))

(defn bray-curtis-dist
  "Same input as adist.  xi is the count of OTU i in sample x."
  [xs ys]
  (assert (= (count xs) (count ys)))
  (- 1
     (/ (* 2
           (apply + (map-indexed (fn [i x]
                                   (let [y (ys i)]
                                     (if (< x y)
                                       x
                                       y)))
                                 xs)))
        (+ (apply + xs)
           (apply + ys)))))

;;;; Data

(def otus [:otu-1 :otu-2 :otu-3])

(def protocol-steps [:extraction-bais
                     :pcr-bias
                     :sequencing-bias
                     :bioinformatics-bais])

;; row i is OTU i, col j is sample j
(def otu-counts
  (r/atom [[1 15] ; OTU 1
           [1  1] ; OTU 2
           [1  4] ; OTU 3
           ]))

;; row i is OTU i, col j is protocol step j
(def protocol-bias
  (r/atom [[ 1 1 1    1]   ; OTU 1
           [ 4 6 0.5  1.5] ; OTU 2
           [15 2 0.25 0.8] ; OTU 3
           ]))

(def num-otus (count @otu-counts))
(def num-samples (count (first @otu-counts)))
(def num-protocol-steps (count protocol-steps))

(def total-bar-height 350)
(def otu-colors ["#004488" "#DDAA33" "#BB5566"])

;; Top left coord for each of the relative abundance bars.
(def rabund-bar-xy-start
  [{:x  96 :y 36.75} ; sample 1 actual
   {:x 216 :y 36.75} ; sample 1 observed
   {:x 384 :y 36.75} ; sample 2 actual
   {:x 504 :y 36.75} ; sample 2 observed
   ])

(def rabund-bar-width 72)
(def svg-width 625)
(def svg-height 550)

;; (def protocol-otu-bias-total
;;   (r/atom nil))

(defn total-bias
  "Returns a vec with the total bias for each OTU.  Deref the atoms first."
  [protocol-bias]
  (vec (map #(reduce * %) protocol-bias)))

(defn observed-counts
  "Returns vec of vecs to match the input data format.  Deref the atoms first."
  [actual-counts protocol-bias]
  (vec (map (fn [bias counts]
              (vec (map #(* bias %) counts)))
            (total-bias protocol-bias)
            actual-counts)))

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
  (defn rabund2 [vec]
    (let [total (apply + vec)]
      (map #(/ % total) vec)))

  (let [;; samples => 2d vec, rows are samples, cols are OTUs
        samples (transpose counts)]
    (transpose (map rabund2 samples))))

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

(defn rabund-bar
  "Make a relative abundance bar for a single sample."
  [sample-idx rect-coords]
  (let [start-xy (get rabund-bar-xy-start sample-idx {:x 0 :y 0})]
    (map-indexed (fn [i coords]
                    (let [x (:x start-xy)
                          ;; offset the y value by the default y start
                          y (+ (:start coords) (:y start-xy))
                          width rabund-bar-width
                          height (:height coords)]
                      ^{:key (str i "-" x "-" y)}
                      [:rect {:x x :y y
                              :width width :height height
                              :fill (otu-colors i)}]))
                  rect-coords)))

;;;; Components

(defn sample-1-labels []
  ;; sample 1 label
  [:g
   [:text {:transform "matrix(1, 0, 0, 1, 192, 490.75)"}
    [:tspan {:x -79.04 :y 12.5
             :font-size 36}
     "Sample 1"]]
   [:g
    ;; sample 1 actual abundance label
    [:text {:transform "matrix(1, 0, 0, 1, 132, 437.25)"}
     [:tspan {:x -19.865 :y -4
              :font-size 13}
      "Actual"]
     [:tspan {:x -34.312 :y 12
              :font-size 13}
      "abundance"]]
    ;; sample 1 observed abundance label
    [:text {:transform "matrix(1, 0, 0, 1, 252, 437.25)"}
     [:tspan {:x -29.986 :y -4
              :font-size 13}
      "Observed"]
     [:tspan {:x -34.312 :y 12
              :font-size 13}
      "abundance"]]]])

(defn sample-2-labels []
  [:g
   ;; sample 2 label
   [:text {:transform "matrix(1, 0, 0, 1, 480, 490.75)"}
    [:tspan {:x -79.04 :y 12.5
             :font-size 36}
     "Sample 2"]]
   [:g
    ;; sample 2 actual abundance label
    [:text {:transform "matrix(1, 0, 0, 1, 420, 437.25)"}
     [:tspan {:x -19.865 :y -4
              :font-size 13}
      "Actual"]
     [:tspan {:x -34.312 :y 12
              :font-size 13}
      "abundance"]]
    ;; sample 2 observed abundance label
    [:text {:transform "matrix(1, 0, 0, 1, 540, 437.25)"}
     [:tspan {:x -29.986 :y -4
              :font-size 13}
      "Observed"]
     [:tspan {:x -34.312 :y 12
              :font-size 13}
      "abundance"]]]])

(defn bar-chart-key []
  [:g
   ;; OTU blue block
   [:rect {:x 5 :y 168
           :width 23 :height 23
           :fill (otu-colors 0)}]
   ;; OTU yellow block
   [:rect {:x 5 :y 204
           :width 23 :height 23
           :fill (otu-colors 1)}]
   ;; OTU red block
   [:rect {:x 5 :y 240
           :width 23 :height 23
           :fill (otu-colors 2)}]
   ;; OTU blue text
   [:text {:transform "matrix(1, 0, 0, 1, 52, 179.5)"}
    [:tspan {:x -20 :y 4} "OTU_1"]]
   ;; OTU yellow text
   [:text {:transform "matrix(1, 0, 0, 1, 52, 216)"}
    [:tspan {:x -20 :y 4} "OTU_2"]]
   ;; OTU red text
   [:text {:transform "matrix(1, 0, 0, 1, 52, 252)"}
    [:tspan {:x -20 :y 4} "OTU_3"]]])

(defn bar-charts
  "Both arguments are atoms."
  [otu-counts otu-bias-per-step]
  (let [actual-rel-abund-by-sample
        (transpose (calc-rel-abund @otu-counts))
        actual-rect-coords
        (map calc-rect-coords actual-rel-abund-by-sample)
        observed-rel-abund-by-sample
        (transpose (calc-rel-abund (observed-counts @otu-counts @protocol-bias)))
        observed-rect-coords
        (map calc-rect-coords observed-rel-abund-by-sample)]
    [:svg {:width svg-width :height svg-height}
     [:g
      (map rabund-bar (range 4) (interleave actual-rect-coords observed-rect-coords))]
     [sample-1-labels]
     [sample-2-labels]
     [bar-chart-key]]))

(defn bias-table-header []
  [:thead
   [:tr
    [:th "OTU ID"]
    [:th "Extraction bias"]
    [:th "PCR bias"]
    [:th "Sequencing bias"]
    [:th "Bioinformatics bias"]
    [:th "Total bias"]]])

#_(defn check-user-bias-input [event]
  (let [val (js/parseFloat (event-val event))]
    (if (or (js/isNaN val) (<= val 0))
      1
      val)))

#_(defn check-user-count-input [event]
  (let [val (js/parseFloat (event-val event))]
    (if (or (js/isNaN val) (< val 1))
      1
      (Math/round val))))

(defn check-user-input [event]
  (let [val (js/parseFloat (event-val event))]
    (if (or (js/isNaN val) (< val 1))
      0
      val)))

(def bias-slider-vals (into [0]
                            (flatten (vector (map #(/ 1 %)
                                                  (reverse (range 1 21)))
                                             (range 2 21)))))
(def bias-slider-min 0)
(def bias-slider-max (count bias-slider-vals))
(def input-focus-bg-color "#ecd6e6")
(def input-blur-bg-color "")

(defn on-focus-change-bg [event]
  (set! (.. event -target -style -background)
        input-focus-bg-color))

(defn bias-table-body
  "Assumes all OTUs have equal number of steps."
  []
  [:tbody
   (doall
    (for [otu-idx (range num-otus)]
      ^{:key (str "otu-" otu-idx)}
      [:tr {:style {:color (otu-colors otu-idx)}}
       [:td (str "OTU_" (inc otu-idx))]
       (doall
        (for [step-idx (range num-protocol-steps)]
          ^{:key (str "step-" step-idx)}
          [:td
           [:input {:type "text"
                    :defaultValue (get-in @protocol-bias
                                          [otu-idx step-idx])
                    :on-focus on-focus-change-bg
                    :on-key-down (fn [event]
                                   (if (= (.-key event) "Enter")
                                     (.. event -target blur)))
                    :on-blur (fn [event]
                               (let [val (check-user-input event)]
                                 (set! (.. event -target -style -background)
                                       input-blur-bg-color)
                                 ;; also set the input value to this thing too...
                                 (set! (.. event -target -value) val)
                                 (swap! protocol-bias
                                        assoc-in
                                        [otu-idx step-idx]
                                        val)))}]]))
       [:td (gstring/format "%.2fx" (reduce * (get @protocol-bias otu-idx)))]]))])

(defn bias-table []
  [:table
   [bias-table-header]
   [bias-table-body]])

(defn count-table-header []
  [:thead
   [:tr
    [:th "OTU ID"]
    [:th "Sample 1"]
    [:th "sample 2"]]])

(defn count-table-body
  "editable: flag to specify whether table data cells should be inputs
  or not. "
  [counts editable]
  [:tbody
   (doall
    (for [otu-idx (range num-otus)]
      ^{:key (str "otu-" otu-idx)}
      [:tr {:style {:color (otu-colors otu-idx)}}
       [:td (str "OTU_" (inc otu-idx))]
       (doall
        (for [sample-idx (range num-samples)]
          ^{:key (str "otu-sample-" otu-idx "-" sample-idx)}
          [:td
           (if (= :editable editable)
             [:input {:type "text"
                      :defaultValue (get-in @counts [otu-idx sample-idx])
                      :on-focus on-focus-change-bg
                      :on-key-down (fn [event]
                                     (if (= (.-key event) "Enter")
                                       (.. event -target blur)))
                      :on-blur (fn [event]
                                 (let [val (Math/round (check-user-input event))]
                                   (set! (.. event -target -style -background)
                                         input-blur-bg-color)
                                   (set! (.. event -target -value) val)
                                   (swap! counts
                                          assoc-in
                                          [otu-idx sample-idx]
                                          val)))}]
             (Math/round (get-in @counts [otu-idx sample-idx])))]))]))])



(defn actual-count-table []
  [:table
   [count-table-header]
   [count-table-body otu-counts :editable]])

(defn observed-count-table []
  [:table
   [count-table-header]
   [count-table-body (r/atom
                      (observed-counts @otu-counts
                                       @protocol-bias))
    :not-editable]])

(defn beta-table [actual-counts observed-counts]
  (let [sample-actual-counts (transpose actual-counts)
        sample-observed-counts (transpose observed-counts)
        s1-actual (first sample-actual-counts)
        s2-actual (last sample-actual-counts)
        s1-observed (first sample-observed-counts)
        s2-observed (last sample-observed-counts)
        actual-euclidean (euclidean-dist s1-actual s2-actual)
        actual-bray-curtis (bray-curtis-dist s1-actual s2-actual)
        actual-aitchison (aitchison-dist s1-actual s2-actual)
        observed-euclidean (euclidean-dist s1-observed s2-observed)
        observed-bray-curtis (bray-curtis-dist s1-observed s2-observed)
        observed-aitchison (aitchison-dist s1-observed s2-observed)]
    [:table
     [:thead
      [:tr
       [:th "S1 vs S2"]
       [:th "Euclidean"]
       [:th "Bray-Curtis"]
       [:th "Aitchison"]]]
     [:tbody
      [:tr
       [:td "Actual"]
       [:td (gstring/format "%.2f" actual-euclidean)]
       [:td (gstring/format "%.2f" actual-bray-curtis)]
       [:td (gstring/format "%.2f" actual-aitchison)]]
      [:tr
       [:td "Observed"]
       [:td (gstring/format "%.2f" observed-euclidean)]
       [:td (gstring/format "%.2f" observed-bray-curtis)]
       [:td (gstring/format "%.2f" observed-aitchison)]]]]))

(defn app-scaffold []
  [:div
   [:div {:class "row"}
    [:div {:class "five columns"}
     [:h3 "Actual counts"]
     [:div {:id "mgs-bias-otu-table-actual-counts"}]]
    [:div {:class "seven columns"}
     [:h3 "Protocol biases"]
     [:div {:id "mgs-bias-protocol-steps"}]]]
   [:div {:class "row"}
    [:div {:class "five columns"}
     [:h3 "Observed counts"]
     [:div {:id "mgs-bias-otu-table-observed-counts"}]
     [:h3 "Sample distance"]
     [beta-table @otu-counts (observed-counts @otu-counts
                                              @protocol-bias)]]
    [:div {:class "seven columns"}
     [:div
      [:h3 "Composition charts"]
      [:div {:id "mgs-bias-composition-charts"}]]]]])

;;;; Rendering

(defn render-app []
  (r/render [app-scaffold]
            (.getElementById js/document html-id-app))
  (r/render [bias-table protocol-bias]
            (.getElementById js/document html-id-bias-table))
  (r/render [actual-count-table]
            (.getElementById js/document html-id-otu-table-actual-counts))
  (r/render [observed-count-table]
            (.getElementById js/document
                             html-id-otu-table-observed-counts))
  (r/render [bar-charts otu-counts protocol-bias]
            (.getElementById js/document html-id-composition-charts)))

(render-app)
