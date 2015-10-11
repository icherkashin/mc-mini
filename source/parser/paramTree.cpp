#include <iostream>
#include <string>
#include <map>

#include "parser/parserException.h"
#include "parser/paramTree.h"

using namespace std;

ParamNode::~ParamNode () {
  for (map<string, ParamNode *>::iterator it = children.begin(); it != children.end(); ++it) {
    delete (it->second);
  }
}

ParamTree::ParamTree () {
  rootNode = focusNode = new ParamNode;
  focusNode->parent = NULL;
  cerr << "Initialized parameter tree." << endl;
}

ParamTree::~ParamTree () {
  delete rootNode;
}

void ParamTree::moveUp (string key) {
  if (focusNode->children.count (key)) {
    map<string, ParamNode *>::iterator iter = focusNode->children.find (key);

    focusNode = iter->second;
  } else {
    throw child_not_found;
  } 
}

void ParamTree::moveDown () {
  if (focusNode->parent == NULL) {
    throw -1;
  } else {
    focusNode = focusNode->parent;
  }
}

string ParamTree::getParam (string key) {
  if (focusNode->params.count (key)) {
    return focusNode->params[key];
  } else {
    throw param_not_found;
  }
}

void ParamTree::addNode (string key) {
  if (focusNode->children.count (key)) {
    throw -1;
  } else {
    ParamNode * newNode = new ParamNode;
    newNode->parent = focusNode;
    focusNode->children.insert(pair<string, ParamNode *> (key, newNode));
  }
}

void ParamTree::delNode (string key) {
  if (focusNode->children.count (key)) {
    ParamNode * match = focusNode->children.find (key)->second;
    focusNode->children.erase (key);
    delete match;
  } else {
    throw child_not_found;
  }
}

void ParamTree::addParam (string key, string param) {
  focusNode->params.insert(pair<string, string> (key, param));
  cerr << "Wrote " << key << " as " << param << endl;
}

void ParamTree::delParam (string key) {
  if (focusNode->params.count (key)) {
    focusNode->params.erase (key);
  } else {
    throw param_not_found;
  }
}
