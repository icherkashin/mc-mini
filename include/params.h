#pragma once

#include "params/param_node.h"

class Params {
public:
  Params(ParamNode *_tree_root);

  template <typename T>
  const T get_param(std::string _key);

  template <typename T>
  const T query_param(std::string _key, const T _default_value) {
    try {
      return get_param<T>(_key);
    } catch(std::invalid_argument e) {
      return _default_value;
    }
  }

  void push(std::string _section_name);
  void try_push(std::string _section_name);

  void pop();

private:
  ParamNode *__root;
  ParamNode *__cursor;
};
