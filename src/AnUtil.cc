// Copyright 2023 nicolò salimbeni andrea de vita
#include "../include/AnUtil.h"

void AnUtil::ProgressBarr(float progress, int present_bar, int total_bars) {
  int barWidth = 70;
  std::cout << present_bar << "/" << total_bars << " [";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos)
      std::cout << "=";
    else if (i == pos)
      std::cout << ">";
    else
      std::cout << " ";
  }
  std::cout << "] " << int(progress * 100.0) << " %\r";
  std::cout.flush();
}
