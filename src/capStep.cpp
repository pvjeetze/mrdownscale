#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

// One cap step of toolAllocateByCrossEntropy, in place (flows and the matrices of capFactors are changed): for each
// cap block and each of its groups, the flows with sign 1 (sum plus) are scaled by a factor a and the flows with
// sign -1 (sum minus) by 1 / a, so that
// plus * a - minus / a = cap, in the grid cells where the group has flows and is at or above its cap or was scaled
// down before. Earlier scaling is reversed up to the cap; the accumulated factor stays at least minFactor.
// flows: active grid cells x flows; cells: row of each active grid cell in the matrices of caps and capFactors;
// groups: per block a list of groups with the cap column and the flow columns with sign 1 and -1 (all 1-based).
// [[Rcpp::export]]
void capStep(NumericMatrix flows, IntegerVector cells, List caps, List capFactors, List groups,
             NumericVector minFactor) {
  for (int b = 0; b < caps.size(); b++) {
    NumericMatrix cap = caps[b];
    NumericMatrix capFactor = capFactors[b];
    List blockGroups = groups[b];
    for (int k = 0; k < blockGroups.size(); k++) {
      List group = blockGroups[k];
      int column = as<int>(group["column"]) - 1;
      IntegerVector plusFlows = group["plus"];
      IntegerVector minusFlows = group["minus"];
      int n = flows.nrow();
      // sums of the flows with sign 1 and -1 per grid cell, column by column and in long double, as rowSums in R
      std::vector<long double> plusSum(n, 0.0L);
      std::vector<long double> minusSum(n, 0.0L);
      for (int j : plusFlows) {
        for (int i = 0; i < n; i++) plusSum[i] += flows(i, j - 1);
      }
      for (int j : minusFlows) {
        for (int i = 0; i < n; i++) minusSum[i] += flows(i, j - 1);
      }

      std::vector<double> factor(n, 1.0);
      bool scaled = false;
      for (int i = 0; i < n; i++) {
        int row = cells[i] - 1;
        double plus = plusSum[i];
        double minus = minusSum[i];
        double capValue = cap(row, column);
        double factorSoFar = capFactor(row, column);
        // only grid cells at or above the cap, or scaled down before, can get a factor other than 1
        if (!(plus + minus > 0 && (plus - minus >= capValue || factorSoFar < 1))) continue;

        // root of plus * a^2 - cap * a - minus = 0, in the form that is exact for cap < 0
        double root = std::sqrt(capValue * capValue + 4 * plus * minus);
        double a = R_PosInf;
        if (capValue < 0) {
          a = 2 * minus / (root - capValue);
        } else if (plus > 0) {
          a = (capValue + root) / (2 * plus);
        }
        // capFactor can reach 0 for hard caps (minFactor 0) where the caps cannot hold
        double lower = minFactor[b] > 0 ? minFactor[b] / factorSoFar : 0;
        a = std::max(std::min(a, 1 / factorSoFar), lower);
        // groups of hard caps with a zero capFactor whose cap does not bind stay as they are
        if (std::isinf(a) || a == 1) continue;
        factor[i] = a;
        capFactor(row, column) *= a;
        scaled = true;
      }
      if (!scaled) continue;

      // a zero factor occurs only where the flows with sign -1 are zero; a factor of 1 leaves a flow unchanged
      for (int j : plusFlows) {
        for (int i = 0; i < n; i++) flows(i, j - 1) *= factor[i];
      }
      for (int j : minusFlows) {
        for (int i = 0; i < n; i++) flows(i, j - 1) *= factor[i] == 0 ? 1 : 1 / factor[i];
      }
    }
  }
}
