#ifndef CODE_SCRIPT_TEMPLATE_H
#define CODE_SCRIPT_TEMPLATE_H

#include "script.h"

namespace Code {

class ScriptTemplate {
public:
  constexpr ScriptTemplate(const char * name, const char * value) : m_name(name), m_value(value) {}
  // INSERT FROM OUTPUT.TXT HERE

static const ScriptTemplate * searcher();

static const ScriptTemplate * importer();

static const ScriptTemplate * AlternatingCurrent();

static const ScriptTemplate * BFieldOfCoilOfWireSolenoid();

static const ScriptTemplate * Diffraction();

static const ScriptTemplate * DopplerEffect();

static const ScriptTemplate * ElectricField();

static const ScriptTemplate * ElectricPotential();

static const ScriptTemplate * ElectromagneticInduction();

static const ScriptTemplate * ElectromagneticSpectrum();

static const ScriptTemplate * Electronvolt();

static const ScriptTemplate * FaradayLaw();

static const ScriptTemplate * Generator();

static const ScriptTemplate * GravitationalAcceleration();

static const ScriptTemplate * HallEffect();

static const ScriptTemplate * HygensPrinciple();

static const ScriptTemplate * InterferenceWaves();

static const ScriptTemplate * LenzLaw();

static const ScriptTemplate * LightWaves();

static const ScriptTemplate * LorentzForce();

static const ScriptTemplate * LrOrRlCircuits();

static const ScriptTemplate * MagneticFields();

static const ScriptTemplate * MagneticFlux();

static const ScriptTemplate * MassSpectrometry();

static const ScriptTemplate * MotionalEmf();

static const ScriptTemplate * NumberOfMaxima();

static const ScriptTemplate * OerstedsExperiment();

static const ScriptTemplate * RefractionOfLight();

static const ScriptTemplate * RlcCircuit();

static const ScriptTemplate * SelfInduction();

static const ScriptTemplate * SlitsWaves();

static const ScriptTemplate * SoundWaves();

static const ScriptTemplate * StandingWaves();

static const ScriptTemplate * Superposition();

static const ScriptTemplate * Torque();

static const ScriptTemplate * TypesOfWaves();

static const ScriptTemplate * WaveEquation();

static const ScriptTemplate * WavesFeatures();

static const ScriptTemplate * WavesInRippleTank();

static const ScriptTemplate * WavesGrating();

static const ScriptTemplate * Youngsdoubleslit();

  // ---------------------------------------------------------------------------------------------

  static const ScriptTemplate * Empty();
  const char * name() const { return m_name; }
  const char * content() const { return m_value + Script::StatusSize(); }
  const char * value() const { return m_value; }
private:
  const char * m_name;
  const char * m_value; // holds the 'importation status' and 'current importation status' flags concatenated with the script content
};

}

#endif
