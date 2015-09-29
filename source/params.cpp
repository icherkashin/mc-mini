#include <iostream>

#include "params.h"

Params::Params(ParamNode *_tree_root) :
    __root(_tree_root), __cursor(_tree_root) {}

/**
 * This preprocessor macro allows us to register a parameter
 * conversion method by providing the type and the function used
 * to convert. The argument `param_t` must be a type, and
 * `conversion_func` must be of the type `const param_t (*)(std::string)`
 */
#define REGISTER_PARAM_CONVERSION(param_t, conversion_func) \
  template <> \
  const param_t Params::get_param<>(std::string _key) { \
    return conversion_func(__cursor->get_data(_key)); \
  }

REGISTER_PARAM_CONVERSION(int, std::stoi);
REGISTER_PARAM_CONVERSION(long, std::stol);
REGISTER_PARAM_CONVERSION(long long, std::stoll);
REGISTER_PARAM_CONVERSION(unsigned long long, std::stoull);
REGISTER_PARAM_CONVERSION(float, std::stof);
REGISTER_PARAM_CONVERSION(double, std::stod);
REGISTER_PARAM_CONVERSION(long double, std::stold);

const std::string string_id(std::string _rvalue) { return _rvalue; }
REGISTER_PARAM_CONVERSION(std::string, string_id)


void Params::push(std::string _section_key) {
  ParamNode *section = __cursor->get_child(_section_key);

  // Make sure that the section actually exists
  if (section == NULL) throw("error");

  // Move the cursor to the section
  __cursor = section;
}

void Params::try_push(std::string _section_key) {
  try {
    // First try pushing normally
    push(_section_key);
  } catch(...) {
    // If pushing failed, create a new dummy node for the section.
    ParamNode *dummy = new ParamNode(__cursor, true);
    __cursor = dummy;
  }
}

void Params::pop() {
  // Ensure that we're not trying to pop from the root
  ParamNode *parent = __cursor->get_parent();

  // Throw an error if we're trying to pop from the root.
  if (parent == NULL) throw("error");

  // Destroy the current node if it's a dummy node.
  if (__cursor->is_dummy()) __cursor->~ParamNode();

  // Move the cursor to the parent node.
  __cursor = parent;
}
