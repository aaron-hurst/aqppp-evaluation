#pragma once
#include <cstdio>
#include <iostream>
#include <vector>

#include "common_content.h"
#include "tool.h"

namespace aqppp {
class Sampling {
 private:
  const double SAMPLE_RATE_;
  const double CI_INDEX_;

 public:
  Sampling(double sample_rate, double ci_index);

  // Given full sample and the demand of user query, cpt the estimate sum result
  // of all data and corresponding confidence interval.
  const std::pair<double, double> Sum(
      const std::vector<std::vector<double>> &sample,
      const std::vector<Condition> &demands) const;
};
}  // namespace aqppp
