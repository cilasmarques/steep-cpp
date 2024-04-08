#include "utils.h"

bool approximatelyEqual(double a, double b)
{
  return fabs(a - b) <= ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * EPS);
}

bool essentiallyEqual(double a, double b)
{
  return fabs(a - b) <= ((fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * EPS);
}

bool definitelyGreaterThan(double a, double b)
{
  return (a - b) > ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * EPS);
}

bool definitelyLessThan(double a, double b)
{
  return (b - a) > ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * EPS);
}

void printVector2x2(const vector<vector<double>> &matrix)
{
  for (const auto &linha : matrix)
  {
    for (const auto &elemento : linha)
    {
      std::cout << elemento << " ";
    }
    std::cout << std::endl;
  }
}
