;;;; Open this link to see the tests: http://localhost:9500/figwheel-extra-main/auto-testing

(ns mgs-bias.core-test
  (:require [cljs.test :refer-macros [deftest is testing run-tests]]
            [mgs-bias.core :as core]))

(deftest say-hi-test
  (is (= "Hi from core.cljs!" (core/say-hi))))

(run-tests)
