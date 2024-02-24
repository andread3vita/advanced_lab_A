// Copyright 2023 nicolò salimbeni andrea de vita
#include "../include/AnUtil.h"
#include <TVectorDfwd.h>

#include <array>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

std::string AnUtil::exec_in_terminal(std::string cmd)
{
  // this function executes in the terminal the command in cmd and saves
  // the ouput in a single string.
  //
  std::array<char, 128>                    buffer;
  std::string                              result;
  std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd.c_str(), "r"), pclose);
  if (!pipe)
  {
    throw std::runtime_error("popen() failed!");
  }
  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
  {
    result += buffer.data();
  }
  return result;
}

void AnUtil::ProgressBar(float progress, std::string title, int present_bar = 1, int total_bars = 1, bool reduce_prints = true)
{
  /*
   arguments explenation: progress as a fraction, ex. 0.3 will be 30% ecc..
                          present bar is an integer
                          total bar is the total number of progress bars in the software

   This function print on screen a progress bar like:

   [1/3] title [==========               ]38%

   Pay attention to the length of the bar.
   This does not fit the current window but uses a fixed length
   */

  // to make is faster print only if the % is an integer
  if (100 * progress - (int)(progress * 100) != 0 && reduce_prints)
  {
    return;
  }

  int barWidth = 70;
  std::cout << "\rTotal (" << present_bar << "/" << total_bars << ") ";
  std::cout << std::setw(30) << title << " [";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i)
  {
    if (i < pos)
      std::cout << "=";
    else if (i == pos)
      std::cout << ">";
    else
      std::cout << " ";
  }
  std::cout << "] " << int(progress * 100.0) << " %";
  if (progress == 1)
    std::cout << std::endl;
  std::cout.flush();
}

TVectorD AnUtil::LoadVector(const std::string &filename, size_t columnNumber)
{
  std::ifstream file(filename);
  if (!file)
  {
    std::cerr << "Failed to open the file." << std::endl;
    return {};
  }

  std::vector<std::vector<double>> columns; // Vector to store columns
  std::string                      line;
  size_t                           numColumns = 0;
  while (std::getline(file, line))
  {
    if (line.empty() || line[0] == '#')
    {
      continue; // Skip empty lines and lines starting with "#"
    }

    std::istringstream  iss(line);
    std::vector<double> row;
    double              value;
    while (iss >> value)
    {
      row.push_back(value);
    }

    if (numColumns == 0)
    {
      numColumns = row.size();
    }
    else if (row.size() != numColumns)
    {
      std::cerr << "Error: Inconsistent number of columns." << std::endl;
      file.close();
      return {};
    }

    columns.push_back(row);
  }

  file.close();

  if (columnNumber - 1 >= numColumns)
  {
    std::cerr << "Error: Column number out of range." << std::endl;
    return {};
  }

  std::vector<double> columnVector;
  for (const auto &row : columns)
  {
    columnVector.push_back(row[columnNumber - 1]);
  }

  TVectorD v_out(columnVector.size());
  for (int i = 0; i < columnVector.size(); i++)
  {
    v_out[i] = columnVector[i];
  }

  return v_out;
}
