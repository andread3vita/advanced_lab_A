#include "./../include/FolderManager.h"
#include <cstdlib>
#include <dirent.h>
#include <iostream>
#include <string>
FolderManager::FolderManager(){};

FolderManager::FolderManager(std::string path)
{
  LoadFolder(path);
}

FolderManager::~FolderManager(){};

std::string FolderManager::GetObjectPath(std::string pattern) const
{
  /*
   This function return the absolute path of the file that contains in the
   name the string pattern (it can be the all name or only a substring).
   But before doing this it checks if there are more than one file whose
   satisfying the request, in that case it returns an errror
   */

  // check if there are more than one file name containing the string pattern
  // in that case return an error
  CheckIfLoaded();
  int         n_files = 0;
  std::string targetFile;
  for (std::string fileName : objects)
  {
    if (fileName.find(pattern) != std::string::npos)
    {
      targetFile = fileName;
      n_files++;
    }
  }
  if (n_files > 1)
  {
    std::cerr << "Error: trying to get the path of a file containing the substring: " + pattern << std::endl;
    std::cerr << "\t---> There a " << n_files << " satisfying this condition, just one needed, abort." << std::endl;
    exit(EXIT_FAILURE);
  }
  else if (n_files == 0)
  {
    std::cerr << "Error: trying to get the path of a file containing the substring: " + pattern << std::endl;
    std::cerr << "\t---> There are " << n_files << " files satisfying this condition, at least one is needed, abort." << std::endl;
    exit(EXIT_FAILURE);
  }

  // if the code gets here it means that there is just one file and it can be returned;
  return pathToFolder + targetFile;
}

std::vector<std::string> FolderManager::GetListOfObjectsPath(std::string pattern) const
{
  // This function return a vector with all the paths of the files whose contain in the name the substring pattern
  // If there are 0 files satisfying the condition it stops the macro.
  CheckIfLoaded();
  std::vector<std::string> targetFiles;

  if (!CheckIfPresent(pattern))
  {
    std::cerr << "ERROR: Trying to get all the files in " + pathToFolder + " containing the substring " + pattern << std::endl;
    std::cerr << "\t---> There are no files with \"" + pattern + "\" in the name, abort." << std::endl;
    exit(EXIT_FAILURE);
  }

  for (std::string fileName : objects)
  {
    if (fileName.find(pattern) != std::string::npos)
    {
      targetFiles.push_back(pathToFolder + fileName);
    }
  }
  return targetFiles;
}

std::vector<std::string> FolderManager::GetAllObjectsPath() const
{
  // This function return all the paths to the objects in the folder in a vector
  CheckIfLoaded();
  std::vector<std::string> allPaths;
  for (int i = 0; i < objects.size(); i++)
  {
    allPaths.push_back(pathToFolder + objects[i]);
  }
  return allPaths;
}

bool FolderManager::CheckIfPresent(std::string pattern) const
{
  // This function checks if there is at least one file with the substring pattern in its name
  CheckIfLoaded();
  bool isPresent = false;
  for (std::string fileName : objects)
  {
    if (fileName.find(pattern) != std::string::npos)
    {
      return true;
    }
  }
  return false;
}

void FolderManager::Print() const
{
  // This function prints on screen the content of the folder
  CheckIfLoaded();

  if (objects.size() == 0 && loaded)
  {
    std::cout << "\nPrint the content of the folder " + pathToFolder << std::endl;
    std::cout << "The folder is empty!\n" << std::endl;
    return;
  }
  else if (!loaded)
  {
    std::cout << "Trying to dump the content of a folder on screen, but the folder is not loaded yet, abort." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout << "\nPrint the content of the folder " + pathToFolder << std::endl;
  for (std::string fileName : objects)
  {
    std::cout << "\t" + fileName << std::endl;
  }
  std::cout << std::endl;
  return;
}

void FolderManager::LoadFolder(std::string path)
{
  // This function load ONLY the files name, not the folders name,
  // contained in the folder in path "path" inside the object vector.
  loaded = true;

  pathToFolder = path;
  if (path[path.size() - 1] != '/')
  {
    path.append("/");
  }
  objects.clear();

  DIR *dir = opendir(path.c_str());

  if (dir)
  {
    struct dirent *entry;

    // Read directory entries
    while ((entry = readdir(dir)) != nullptr)
    {
      if (entry->d_type == DT_REG)
      { // Check if it's a regular file
        objects.push_back(entry->d_name);
      }
    }

    closedir(dir);
  }
  else
  {
    std::cerr << "ERROR opening directory: " << path << ", abort." << std::endl;
    exit(EXIT_FAILURE);
  }

  return;
}

std::vector<std::string> FolderManager::GetListOfObjectsPathComposed(std::vector<std::string> patters) const
{
  // This function return a vector with all the paths of the files whose contain in the name the substring pattern
  // If there are 0 files satisfying the condition it stops the macro.
  CheckIfLoaded();
  std::vector<std::string> targetFiles;

  for (std::string fileName : objects)
  {
    bool isOk = true;
    for (std::string pattern : patters)
    {
      if (fileName.find(pattern) == std::string::npos)
      {
        isOk = false;
      }
    }
    if (isOk)
    {
      targetFiles.push_back(pathToFolder + fileName);
    }
  }

  if (targetFiles.size() == 0)
  {
    std::cerr << "ERROR: Trying to get files in " + pathToFolder + " containing the following substrings:" << std::endl;
    for (std::string substring : patters)
    {
      std::cerr << "\t" + substring << std::endl;
    }
    std::cerr << "\tBut there aren't files satisfying these conditions" << std::endl;
    exit(EXIT_FAILURE);
  }
  return targetFiles;
}

std::string FolderManager::GetObjectPathComposed(std::vector<std::string> patters) const
{
  // This function is the same of GetObjectPath, but instead of just one string here
  // I look for a set of substrings. This is usefull because sometimes you have to deal with
  // files with very similar name.

  CheckIfLoaded();
  std::vector<std::string> candidate = GetListOfObjectsPathComposed(patters);

  if (candidate.size() > 1)
  {
    std::cerr << "ERROR: Trying to get a file in " + pathToFolder + " containing the following substrings:" << std::endl;
    for (std::string substring : patters)
    {
      std::cerr << "\t" + substring << std::endl;
    }
    std::cerr << "\tBut there are more than one file satisfying this condition, abort." << std::endl;
    exit(EXIT_FAILURE);
  }
  return candidate[0];
}

void FolderManager::CheckIfLoaded() const
{
  // Check if the folder is loaded, if not stops the software.
  // This function is used at the beginning of other functions as a security check.
  if (!loaded)
  {
    std::cerr << "ERROR: trying to call a function of the class FolderManager on an instance where the folder is not loaded yet, abort" << std::endl;
    exit(EXIT_FAILURE);
  }
}
