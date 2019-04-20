(ns mgs-bias.core
  (:require
   [goog.string :as gstring]
   [goog.string.format]
   [goog.dom :as gdom]
   [reagent.core :as r]))

(defn say-hi []
  "Hi from core.cljs!")

(def html-id-bias-table "mgs-bias-protocol-steps")
(def html-id-num-otus "mgs-bias-num-otus")

(def config (atom {:num-otus 3}))

;;;; Data

(def bias-table-data
  (r/atom {:otu-1
           {:actual-count 1
            :extraction-bias 1
            :pcr-bias 1
            :sequencing-bias 1
            :bioinformatics-bias 1}
           :otu-2
           {:actual-count 1
            :extraction-bias 4
            :pcr-bias 6
            :sequencing-bias 0.5
            :bioinformatics-bias 1.5}
           :otu-3
           {:actual-count 1
            :extraction-bias 15
            :pcr-bias 2
            :sequencing-bias 0.25
            :bioinformatics-bias 0.8}}))

;;;; Doing things to data

(def bias-table-columns 9)

(defn total-bias [otu-info]
  (* (:extraction-bias otu-info)
     (:pcr-bias otu-info)
     (:sequencing-bias otu-info)
     (:bioinformatics-bias otu-info)))

(defn observed-count [otu-info]
  (* (:actual-count otu-info)
     (total-bias otu-info)))

(defn change-table-val [otu-id key new-val]
  (swap! bias-table-data
           assoc-in
           [otu-id key]
           new-val))

(defn event-val [event]
  (-> event .-target .-value))


;;;; Components

(defn bias-table-header []
  [:thead
   [:tr
    [:th "Actual count"]
    [:th "Actual composition"]
    [:th "Extraction bias"]
    [:th "PCR bias"]
    [:th "Sequencing bias"]
    [:th "Bioinformatics bias"]
    [:th "Total bias"]
    [:th "Observed count"]
    [:th "Observed composition"]]])

(defn mgs-bias-table-input [table-data otu-id otu-info key]
  [:input {:type "text"
           :value (key otu-info)
           :on-change #(swap! table-data
                              assoc-in
                              [otu-id key]
                              (event-val %))}])

(defn mgs-bias-protocol-table [data]
  [:table
   [bias-table-header]
   [:tbody
    (map (fn [[id info]]
           ^{:key id}
           [:tr
            [:td [mgs-bias-table-input bias-table-data id info :actual-count]]
            [:td "actual-composition"]
            [:td [mgs-bias-table-input bias-table-data id info :extraction-bias]]
            [:td [mgs-bias-table-input bias-table-data id info :pcr-bias]]
            [:td [mgs-bias-table-input bias-table-data id info :sequencing-bias]]
            [:td [mgs-bias-table-input bias-table-data id info :bioinformatics-bias]]
            [:td (gstring/format "%.2f" (total-bias info))]
            [:td (Math/round (observed-count info))]
            [:td "observed-composition"]])
         @bias-table-data)]])

(defn run []
  (r/render [mgs-bias-protocol-table]
            (.getElementById js/document html-id-bias-table)))

(run)
