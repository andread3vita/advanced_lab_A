#ifndef FolderManager_h
#define FolderManager_h

#include <string>
#include <vector>

/*
This class has the purpose of conveniently handling folders in C++,
specifically written libraries could be used but to keep the code as simple
as possible this library was written. It allows you to have the main
functions as: return the global paths of one or more files, the names of the
objects and have the contents of the folder printed on the screen.
*/

class FolderManager
{
  public:
  FolderManager();
  FolderManager(std::string path);
  ~FolderManager();

  void LoadFolder(std::string path); // This function load a folder in the class

  std::string GetObjectPath(std::string pattern) const;                      // This returns the path to an object in that folder
  std::string GetObjectPathComposed(std::vector<std::string> pattern) const; // This returns the path to an object to a folder that contains in the name more than one substring

  std::vector<std::string> GetListOfObjectsPath(std::string pattern) const;                      // This return a vector with all the paths to the files whose name contain the string pattern
  std::vector<std::string> GetListOfObjectsPathComposed(std::vector<std::string> pattern) const; // This return a vector with all the paths to the files whose name constains a set of substrings

  std::vector<std::string> GetAllObjectsPath() const;                 // This returns all the paths to the files in the folder
  bool                     CheckIfPresent(std::string pattern) const; // This function checks if there is at least one file with the pattern substring in the name
  void                     Print() const;                             // Print on the screen what is inside the folder

  private:
  void CheckIfLoaded() const; // this function checks if the folder was loaded, if not stops the software

  bool                     loaded = false;
  std::string              pathToFolder;
  std::vector<std::string> objects;
};

#endif
