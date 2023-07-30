#ifndef CODE_SCRIPT_TEMPLATE_H
#define CODE_SCRIPT_TEMPLATE_H

#include "script.h"

namespace Code {

class ScriptTemplate {
public:
  constexpr ScriptTemplate(const char * name, const char * value) : m_name(name), m_value(value) {}
  static const ScriptTemplate * Empty();
  static const ScriptTemplate * Squares();
  static const ScriptTemplate * Mandelbrot();
  static const ScriptTemplate * Polynomial();
  static const ScriptTemplate * Parabola();
static const ScriptTemplate * record0_();
static const ScriptTemplate * record1_();
static const ScriptTemplate * record2_();
static const ScriptTemplate * record3_();
static const ScriptTemplate * record4_();
static const ScriptTemplate * record5_();
static const ScriptTemplate * record6_();
static const ScriptTemplate * record7_();
static const ScriptTemplate * record8_();
static const ScriptTemplate * record9_();
static const ScriptTemplate * record10_();
static const ScriptTemplate * record11_();
static const ScriptTemplate * record12_();
static const ScriptTemplate * record13_();

  const char * name() const { return m_name; }
  const char * content() const { return m_value + Script::StatusSize(); }
  const char * value() const { return m_value; }
private:
  const char * m_name;
  const char * m_value; // holds the 'importation status' and 'current importation status' flags concatenated with the script content
};

}

#endif
