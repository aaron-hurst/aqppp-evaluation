// #include "stdafx.h"
#include "sampling.h"

namespace aqppp {

Sampling::Sampling(double sample_rate, double ci_index)
    : SAMPLE_RATE_(sample_rate), CI_INDEX_(ci_index) {}

const std::pair<double, double> Sampling::Sum(
    const std::vector<std::vector<double>> &sample,
    const std::vector<Condition> &demands) const {
  assert(!sample.empty());
  const double N_ROWS = sample[0].size();
  assert(N_ROWS > 1);
  double sum = 0.0, sum2 = 0.0;
  for (int rowi = 0; rowi < sample[0].size(); rowi++) {
    bool flag = true;
    for (int i = 0; i < demands.size(); i++) {
      if (DoubleLess(sample[i + 1][rowi], demands[i].lb) ||
          DoubleGreater(sample[i + 1][rowi], demands[i].ub)) {
        flag = false;
        break;
      }
    }
    if (!flag) continue;
    sum += sample[0][rowi];
    sum2 += sample[0][rowi] * sample[0][rowi];
  }
  double variance = 0;
  variance = sum2 / N_ROWS - (sum / N_ROWS) * (sum / N_ROWS);
  double ci = CI_INDEX_ * sqrt(variance / N_ROWS) * N_ROWS / SAMPLE_RATE_;
  // cout << "Sampling: sample rate, sum2 sum2: " << PAR.SAMPLE_RATE <<"
  // "<<sum2<<"
  // "<<sum<< endl;
  return std::pair<double, double>(sum / SAMPLE_RATE_, ci);
}

}  // namespace aqppp
