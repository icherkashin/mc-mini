#include "params/param_node.h"

ParamNode::ParamNode() :
    __parent(NULL) {}

ParamNode::ParamNode(ParamNode *_parent, bool _is_dummy) :
    __parent(_parent), __is_dummy(_is_dummy) {}

ParamNode *ParamNode::get_child(std::string _child_key) {
  return __child_map[_child_key];
}

std::string ParamNode::get_data(std::string _data_key) {
  return __data_map[_data_key];
}

void ParamNode::set_child(std::string _child_key, ParamNode *_child) {
  __child_map[_child_key] = _child;
}

void ParamNode::set_data(std::string _data_key, std::string _data) {
  __data_map[_data_key] = _data;
}

ParamNode *ParamNode::get_parent() {
  return __parent;
}

const bool ParamNode::is_dummy() {
  return __is_dummy;
}
