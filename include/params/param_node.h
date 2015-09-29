#pragma once

#include <string>
#include <map>

class ParamNode {
public:
  ParamNode();
  ParamNode(ParamNode *_parent, bool _is_dummy = false);

  ParamNode *get_child(std::string _child_key);
  std::string get_data(std::string _data_key);

  void set_child(std::string _child_key, ParamNode *_child);
  void set_data(std::string _data_key, std::string _data);

  ParamNode *get_parent();

  const bool is_dummy();

private:
  ParamNode *__parent;

  bool __is_dummy;

  std::map<std::string, ParamNode*> __child_map;
  std::map<std::string, std::string> __data_map;
};
