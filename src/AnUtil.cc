// Copyright 2023 nicolò salimbeni andrea de vita
#include "../include/AnUtil.h"
#include <TVectorDfwd.h>

#include <array>
#include <cstdio>
#include <fstream>
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

void AnUtil::ProgressBarr(float progress, int present_bar, int total_bars)
{
  int barWidth = 70;
  std::cout << present_bar << "/" << total_bars << " [";
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
  std::cout << "] " << int(progress * 100.0) << " %\r";
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
