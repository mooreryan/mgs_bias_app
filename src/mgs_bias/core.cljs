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

(def mgs-bias-data
  (r/atom {:otu-1
           {:otu :otu-1
            :actual-1 1
            :actual-2 15
            :observed-1 1
            :observed-2 15
            :extraction-bias 1
            :pcr-bias 1
            :sequencing-bias 1
            :bioinformatics-bias 1}
           :otu-2
           {:otu :otu-2
            :actual-1 1
            :actual-2 1
            :observed-1 18
            :observed-2 18
            :extraction-bias 4
            :pcr-bias 6
            :sequencing-bias 0.5
            :bioinformatics-bias 1.5}
           :otu-3
           {:otu :otu-3
            :actual-1 1
            :actual-2 4
            :observed-1 6
            :observed-2 24
            :extraction-bias 15
            :pcr-bias 2
            :sequencing-bias 0.25
            :bioinformatics-bias 0.8}}))

;;;; Doing things to data

(defn composition-charts [data]
  (.log js/console (clj->js data))
  ;; Make sure any existing svg has been removed
  (-> js/d3
      (.select (str "#" html-id-composition-charts-svg))
      .remove)
  ;; Append the svg
  (let [svg (-> js/d3
                (.select (str "#" html-id-composition-charts))
                (.append "svg")
                (.attr "id" html-id-composition-charts-svg)
                (.attr "width" 500)
                (.attr "height" 500))
        actual-1 (map (fn [m]
                        (js/parseFloat (:actual-1 m)))
                      data)
        actual-1-rel (map (fn [n]
                            (/ n (reduce + actual-1)))
                          actual-1)
        rects (-> svg
                  (.append "g")
                  (.selectAll "rect")
                  (.data (clj->js actual-1-rel)))
        actual-2 (map (fn [m]
                        (js/parseFloat (:actual-2 m)))
                      data)
        observed-1 (map (fn [m]
                          (js/parseFloat (:observed-1 m)))
                        data)
        observed-2 (map (fn [m]
                          (js/parseFloat (:observed-1 m)))
                        data)]

    (.log js/console (clj->js actual-1-rel))
    (-> rects
        .enter
        (.append "rect")
        (.merge rects)
        (.attr "fill" (fn [d i] (get ["red" "green" "blue"] i)))
        (.attr "width" 40)
        (.attr "height" (fn [d] (* 400 d)))
        (.attr "x" 10)
        (.attr "y" (fn [d i] (.log js/console (* i (* 400 d))) (* i (* 400 d)))))
    (-> rects
        .exit
        .remove)))

(def bias-table-columns 9)

(defn total-bias [otu-info]
  (* (:extraction-bias otu-info)
     (:pcr-bias otu-info)
     (:sequencing-bias otu-info)
     (:bioinformatics-bias otu-info)))

(defn observed-count [otu-info]
  (* (:actual-count otu-info)
     (total-bias otu-info)))

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

(defn table-input [table-data id key]
  (let [change-fn (fn [event]
                    (swap! table-data
                           assoc-in
                           [id key]
                           (event-val event))
                    (doall
                     (let [bias (total-bias (id @table-data))]
                       (for [i [1 2]]
                         (swap! table-data
                                assoc-in
                                [id (keyword (str "observed-" i))]
                                (* bias
                                   (js/parseFloat (get-in @table-data
                                                          [id (keyword (str "actual-" i))])))))))
                    (.log js/console "mgs-bias-data" (clj->js @mgs-bias-data))
                    (composition-charts (vals @mgs-bias-data)))]
    [:input {:type "text"
             :value (get-in @table-data [id key])
             :on-change change-fn}]))

(defn mgs-bias-protocol-table [data]
  [:table
   [bias-table-header]
   [:tbody
    (doall
     (for [otu-id (keys @data)]
       ^{:key otu-id}
       [:tr
        [:td otu-id]
        [:td [table-input data otu-id :extraction-bias]]
        [:td [table-input data otu-id :pcr-bias]]
        [:td [table-input data otu-id :sequencing-bias]]
        [:td [table-input data otu-id :bioinformatics-bias]]
        [:td (gstring/format "%.2f" (total-bias (otu-id @data)))]]))]])

(defn otu-table-actual [data]
  [:table
   [:thead
    [:tr
     [:th "OTU ID"]
     [:th "Sample 1"]
     [:th "sample 2"]]]
   [:tbody
    (doall
     (for [otu-id (keys @data)]
       ^{:key otu-id}
       [:tr
        [:td otu-id]
        [:td [table-input data otu-id :actual-1]]
        [:td [table-input data otu-id :actual-2]]]))]])

(defn otu-table-observed [data]
  [:table
   [:thead
    [:tr
     [:th "OTU ID"]
     [:th "Sample 1"]
     [:th "sample 2"]]]
   [:tbody
    (doall
     (for [otu-id (keys @data)]
       ^{:key otu-id}
       [:tr
        [:td otu-id]
        [:td (:observed-1 (otu-id @data))]
        [:td (:observed-2 (otu-id @data))]]))]])

;;;; Rendering functions

(defn render-mgs-bias-protocol-table [html-id data]
  (r/render [mgs-bias-protocol-table data]
            (.getElementById js/document html-id)))

(defn render-otu-table-actual [html-id data]
  (r/render [otu-table-actual data]
            (.getElementById js/document html-id)))

(defn render-otu-table-observed [html-id data]
  (r/render [otu-table-observed data]
            (.getElementById js/document html-id)))



(render-mgs-bias-protocol-table html-id-bias-table mgs-bias-data)
(render-otu-table-actual html-id-otu-table-actual-counts mgs-bias-data)
(render-otu-table-observed html-id-otu-table-observed-counts mgs-bias-data)

(composition-charts (vals @mgs-bias-data))



;;;; scratch

(defn stack [keys]
  (-> js/d3
      .stack
      (.keys (clj->js keys))
      (.order js/d3.stackOrderNone)
      (.offset js/d3.stackOffsetNone)))

((stack [:actual-1 :observed-1 :actual-2 :observed-2])
 (clj->js (vals @mgs-bias-data)))
