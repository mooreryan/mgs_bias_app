;;;; Open this link to see the tests: http://localhost:9500/figwheel-extra-main/auto-testing

(ns mgs-bias.core-test
  (:require [cljs.test :refer-macros [deftest is testing run-tests]]
            [mgs-bias.core :as core]
            [goog.string :as gstring]))

(deftest say-hi-test
  (is (= "Hi from core.cljs!" (core/say-hi))))

(defn float=
  "Converts to a string rounded to ten-thousandth place then compares
  those."
  [f1 f2]
  (=
   (gstring/format "%.4" f1)
   (gstring/format "%.4" f2)))

(deftest euclidean-dist-test
  (is (float= 0
              (core/euclidean-dist [1 2 3] [1 2 3])))
  (is (float= 37.41657
              (core/euclidean-dist [0 0 0] [10 20 30])))
  (is (float= 33.67492
              (core/euclidean-dist [1 2 3] [10 20 30])))
  (is (float= 33.67492
              (core/euclidean-dist [0 2 3] [10 20 30]))
      "Even one zero value is fine"))

(deftest bray-curtis-dist-test
  (is (float= 0
              (core/euclidean-dist [1 2 3] [1 2 3])))
  (is (float= 1
              (core/euclidean-dist [0 0 0] [10 20 30])))
  (is (float= 0.8181818
              (core/euclidean-dist [1 2 3] [10 20 30])))
  (is (float= 0.8461538
              (core/euclidean-dist [0 2 3] [10 20 30]))
      "Even one zero value is fine"))

(run-tests)
