#include "debug.h"

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