#include "script_template.h"

namespace Code {

// INSERT FROM OUTPUT.TXT HERE

constexpr ScriptTemplate searcherScriptTemplate("searcher.py", "\x01" R"(from importer import * 

vars_names=[AlternatingCurrent,BFieldOfCoilOfWireSolenoid,Diffraction,DopplerEffect,ElectricField,ElectricPotential,ElectromagneticInduction,ElectromagneticSpectrum,Electronvolt,FaradayLaw,Generator,GravitationalAcceleration,HallEffect,HygensPrinciple,InterferenceWaves,LenzLaw,LightWaves,LorentzForce,LrOrRlCircuits,MagneticFields,MagneticFlux,MassSpectrometry,MotionalEmf,NumberOfMaxima,OerstedsExperiment,RefractionOfLight,RlcCircuit,SelfInduction,SlitsWaves,SoundWaves,StandingWaves,Superposition,Torque,TypesOfWaves,WaveEquation,WavesFeatures,WavesInRippleTank,WavesGrating,Youngsdoubleslit]

def s(in_val,characters=20):
    global vars_names
    for text in vars_names:
        find_res=text.find(in_val)
        while find_res!=-1:
            print("IN "+text.split(' ')[0]+'\n'+text[max(0,find_res-characters):min(find_res + characters, len(text))]+'\n')
            find_res=text.find(in_val,find_res+1))");

const ScriptTemplate * ScriptTemplate::searcher() {return &searcherScriptTemplate;}

constexpr ScriptTemplate importerScriptTemplate("importer.py", "\x01" R"(from AlternatingCurrent import * 
from BFieldOfCoilOfWireSolenoid import * 
from Diffraction import * 
from DopplerEffect import * 
from ElectricField import * 
from ElectricPotential import * 
from ElectromagneticInduction import * 
from ElectromagneticSpectrum import * 
from Electronvolt import * 
from FaradayLaw import * 
from Generator import * 
from GravitationalAcceleration import * 
from HallEffect import * 
from HygensPrinciple import * 
from InterferenceWaves import * 
from LenzLaw import * 
from LightWaves import * 
from LorentzForce import * 
from LrOrRlCircuits import * 
from MagneticFields import * 
from MagneticFlux import * 
from MassSpectrometry import * 
from MotionalEmf import * 
from NumberOfMaxima import * 
from OerstedsExperiment import * 
from RefractionOfLight import * 
from RlcCircuit import * 
from SelfInduction import * 
from SlitsWaves import * 
from SoundWaves import * 
from StandingWaves import * 
from Superposition import * 
from Torque import * 
from TypesOfWaves import * 
from WaveEquation import * 
from WavesFeatures import * 
from WavesInRippleTank import * 
from WavesGrating import * 
from Youngsdoubleslit import * )");

const ScriptTemplate * ScriptTemplate::importer() {return &importerScriptTemplate;}

constexpr ScriptTemplate AlternatingCurrentScriptTemplate("AlternatingCurrent.py", "\x01" R"(AlternatingCurrent = """AlternatingCurrent
Alternating Current (AC) is
a type of electric current
that periodically reverses
direction. It is the most
common type of electrical
current used in homes,
businesses, and industrial
settings. The direction of
AC changes approximately 50
or 60 times per second,
depending on the country.
This change in direction
allows AC to be transmitted
over long distances without
significant energy loss.

Io - maximum current
\<F ru> - Phase change

I(t) = Io * sin(2 * pi * f *
t + \<F ru>)
""")");

const ScriptTemplate * ScriptTemplate::AlternatingCurrent() {return &AlternatingCurrentScriptTemplate;}

constexpr ScriptTemplate BFieldOfCoilOfWireSolenoidScriptTemplate("BFieldOfCoilOfWireSolenoid.py", "\x01" R"(BFieldOfCoilOfWireSolenoid = """BFieldOfCoilOfWireSolenoid
aka Solenoid
Magnitnoe pole prohodit ot +
do - vnutri "Coil of wire",
po bokam magnitnoe pole idet
obratno i ogibaet "Coil of
wire". 

Magnetic field in solenoid
directly proportional to the
current

Current in circuit with
solenoid will grow
exponentially, slowly.

L - Inductance [Hanry]
I - Current
t - time
emf: back emf [V]

Emf back = - L * I / t

Inductance:

B = \<miu> * n * I -
Magnetic field in a solenoid
\<F ru> = B * A *
cos(\<theta>) - Magnetic
flux, cos(\<theta>) = 1
\<F ru> = N * \<miu> * n * I
* A = 
= \<miu> * n * \<big i> * A
* N
emf = - \<delta>\<F ru> /
\<delta>t = 
= - (\<miu> * n * A * N) *
\<big i> / t
L = - (\<miu> * A * N ** 2)
/ l

Effect on circuit:
V = I * R
V = emf + IR
""")");

const ScriptTemplate * ScriptTemplate::BFieldOfCoilOfWireSolenoid() {return &BFieldOfCoilOfWireSolenoidScriptTemplate;}

constexpr ScriptTemplate DiffractionScriptTemplate("Diffraction.py", "\x01" R"(Diffraction = """Diffraction
Diffraction: Diffraction is
the phenomenon where waves,
such as light, sound, or
even water waves, spread out
after passing through an
opening or around an
obstacle. It is a result of
the interaction of waves
with the boundaries of the
obstacle or opening. The
degree of spreading, or
diffraction, depends on the
relative sizes of the wave
and the obstacle or opening.

see Hygens_Principle

If an object of width a is
causing the diffraction, the
amount of diffraction can be
determined (qualitatively):

== - ALMOST EQUAL

Very significant
diffraction:
\<lambda> / a == 1 or
\<lambda> > a

Less diffraction:
\<lambda> / a < 1

No diffraction (about 10 **
(-3)):
\<lambda> / a << 1

FROM NOTES:
Single slit diffraction:

\<alpha> - angle between
centre and first minimum
\<lambda> - wavelength
a - slit width

sin(\<alpha>) = 2x / a
x = (n + 0.5) * \<lambda>, n
= 0, x = 1/2 * \<lambda>

sin(\<alpha>) = \<lambda> /
a

y - distance between max and
min on the screen
D - distance between
diffraction object and the
screen

tan(\<alpha>) = y / D

WITH BOOKLET:
Single vs Double slit:
Single:
sin(\<alpha>) = n *
\<lambda> / d
x = L * tan(\<alpha>) -
shows distance between min
and center on the screen

Double:
x = n * \<lambda> * L / d -
shows distance between max
and center on the screen

Diffraction grating:
A diffraction grating is an
optical grating with a
periodic structure that
diffracts light into several
beams traveling in different
directions, i.e., different
diffraction angles. The
directions or diffraction
angles of these beams depend
on the wave (light) incident
angle to the diffraction
grating, the spacing or
distance between adjacent
diffracting elements (e.g.,
parallel slits for a
transmission grating) on the
grating, and the wavelength
of the incident light. The
grating acts as a dispersive
element, meaning it
separates an incident
polychromatic beam into its
constituent wavelength
components at different
angles

\<alpha> - angle between
centre and n-th max
\<lambda> - wavelength
d - separation between slits
n - number of maximum

sin(\<alpha>) = n *
\<lambda> / d

Distance on the screen:
y = L * tan(\<alpha>)
""")");

const ScriptTemplate * ScriptTemplate::Diffraction() {return &DiffractionScriptTemplate;}

constexpr ScriptTemplate DopplerEffectScriptTemplate("DopplerEffect.py", "\x01" R"(DopplerEffect = """DopplerEffect
Doppler effect is the change
in frequency of a wave
detected by any observer
because the wave source and
the observer have different
velocities with respect to
the medium of the wave
propagation.

f' - shifted frequency
v - speed of sound in medium
vs - speed of source
vo - speed of observer
fs - original frequency

Source moves towards
observer: f' = v / (v - vs)
* fs
Source moves against
observer: f' = v / (v + vs)
* fs
Observer moves towards
source: f' = (v + vo) / v *
fs
Observer moves against
source: f' = (v - vo) / v *
fs 
""")");

const ScriptTemplate * ScriptTemplate::DopplerEffect() {return &DopplerEffectScriptTemplate;}

constexpr ScriptTemplate ElectricFieldScriptTemplate("ElectricField.py", "\x01" R"(ElectricField = """ElectricField
Region around a charged
particle in which another
charged particle will
experience a force.


""")");

const ScriptTemplate * ScriptTemplate::ElectricField() {return &ElectricFieldScriptTemplate;}

constexpr ScriptTemplate ElectricPotentialScriptTemplate("ElectricPotential.py", "\x01" R"(ElectricPotential = """ElectricPotential
Electric potential
difference: the change in
electric potential energy of
a positive charge between
two points in a field.
\<delta>V = \<delta>E / q
(Energy per unit of charge)

Determination:
W = \<delta>E = F x s
\<delta>E = ((k * Q * q) / r
** 2) x r
\<delta>E = k * q * Q / r
\<delta>V = \<delta>E / q
V = k * Q / r
Electric potential at a
point a distance r from
charge Q

Definition of electric
potential: The work needed
to bring a positive point
charge into an electric
field from a reference point
where V = 0 (r = \<inf>)
into the field, a distance r
from Q.
V = k * Q / r

Equipotential lines: Lines
where electric potential are
the same

Work needed to move a charge
from one equipotential line
to another (from A to B):
Q - charge, which we are
moving
Vb - electric potential at
the point B
Va - electric potential at
the point A

W = Q * (Vb - Va)
""")");

const ScriptTemplate * ScriptTemplate::ElectricPotential() {return &ElectricPotentialScriptTemplate;}

constexpr ScriptTemplate ElectromagneticInductionScriptTemplate("ElectromagneticInduction.py", "\x01" R"(ElectromagneticInduction = """ElectromagneticInduction
Electromagnetic induction is
the creation of an
electro-motive force (EMF)
by way of a moving magnetic
field around an electric
conductor and, conversely,
the creation of current by
moving an electric conductor
through a static magnetic
field. Electromagnetic
interference (EMI) is also
known as electric current
and electromagnetic
induction and may also be
called magnetic induction,
as the principle remains the
same whether the process is
carried out through
electromagnet or static
magnet.

emf = - \<delta>\<F ru> /
\<delta>t (For a solenoid,
the flux through a single
winding must be multiplied
by the number of windings.)
""")");

const ScriptTemplate * ScriptTemplate::ElectromagneticInduction() {return &ElectromagneticInductionScriptTemplate;}

constexpr ScriptTemplate ElectromagneticSpectrumScriptTemplate("ElectromagneticSpectrum.py", "\x01" R"(ElectromagneticSpectrum = """ElectromagneticSpectrum
- Gamma rays (X - rays)
Highest frequency, emitted
when nuclei make quantum
jumps (e.g. radioactive
decay), highly penetrating,
uses: Induce mutation in
genetic experiments,
sanitize equipment, gamma
ray astronomy. (10 ** 20 ->
10 ** 26 Hz)

- X-rays
Produced when high energy
electrons collide with heavy
metal target, transferring
all kinetic energy to X-ray
photons, uses: x-ray
machine, CT scans, X-ray
astronomy (10 ** 17 -> 10 **
20 Hz)

- Ultraviolet(UV)
Producing: passing electric
current through gas, uses:
body uses UV to produce
vitamin D, UV spectroscopy
in pharmaceutical studies.
(10 ** 15 -> 10 ** 17)

- Visible light
Warm bodies emit visible
light (Black body
radiation), uses: vision in
most living creatures,
Chlorophyll absorbs for
photosynthesis. (4 * 10 **
14 -> 10 ** 15 Hz)

- Infrared
Warm bodies emit IR, uses:
some animals use for
detecting prey, infrared
detectors (security/army) as
they can escape dust. (10 **
13 -> 4 * 10 ** 14)

- Microwaves
Produced by warm bodies,
uses: Transmit TV and
telephone signals,
microwaves to heat food. (10
** 10 -> 10 ** 13 Hz)

- Radio waves
Usually produced by changing
electric fields or changing
magnetic fields, uses: radio
and IV signals. (10 ** 3 ->
10 ** 10 Hz)


""")");

const ScriptTemplate * ScriptTemplate::ElectromagneticSpectrum() {return &ElectromagneticSpectrumScriptTemplate;}

constexpr ScriptTemplate ElectronvoltScriptTemplate("Electronvolt.py", "\x01" R"(Electronvolt = """Electronvolt
Definition: The work done to
accelerate one electron
through a potential
difference of 1 Volt.

W = Q * U
1 eV = 1.6* 10 ** (-19) J

Matter:
E = m * c ** 2

Photon: 
h - plank constant
f - frequency

E = h * f

Photon: Carrier of
electromagnetic wave, acts
like a wave and a particle

""")");

const ScriptTemplate * ScriptTemplate::Electronvolt() {return &ElectronvoltScriptTemplate;}

constexpr ScriptTemplate FaradayLawScriptTemplate("FaradayLaw.py", "\x01" R"(FaradayLaw = """FaradayLaw
Faraday’s law of induction
is a relationship expressing
that a changing magnetic
field induces a voltage in a
circuit.

The electromotive force
or _EMF_ refers to the
potential difference across
the _unloaded_ loop (i.e.
when the resistance in the
circuit is high). In
practice it is often
sufficient to think of EMF
as voltage since both
voltage and EMF are measured
using the same unit,
the volt.


""")");

const ScriptTemplate * ScriptTemplate::FaradayLaw() {return &FaradayLawScriptTemplate;}

constexpr ScriptTemplate GeneratorScriptTemplate("Generator.py", "\x01" R"(Generator = """Generator
A generator works by
converting mechanical energy
into electrical energy. This
conversion process is based
on the principles of
electromagnetic induction,
which were discovered by
physicist Michael Faraday in
1831.
The basic working principle
of a generator involves the
movement of a conductor,
such as a wire, in a
magnetic field. When a coil
of wire (also referred to as
the armature) is rotated in
a magnetic field, an
electric current is induced
in the wire.
The key components of a
generator typically include
a magnet (to create the
magnetic field), an armature
(a coil of wire), and a
commutator. These components
work together to convert
mechanical energy, often
derived from a prime mover
like a diesel engine, water
turbine, or steam turbine,
into electrical energy.
Here's a simplified
explanation of how it works:

1. A magnetic field is
established by the rotating
magnet.
2. The armature (coil of
wire) is placed in this
magnetic field.
3. As the armature rotates,
the magnetic field passes
through the armature,
causing an electromotive
force (EMF) to be induced in
the armature. This is due to
Faraday's law of
electromagnetic induction.
4. The induced EMF causes
electrons to move in the
armature, creating an
electric current.
5. The electric current is
collected and used to power
electrical devices or stored
in batteries.

The output of a generator is
influenced by several
parameters, including the
strength of the magnetic
field, the speed of the
mechanical input (often
represented as the speed of
rotation of the armature),
the number of turns in the
coil, and the coil’s area.
Adjusting these parameters
appropriately can influence
the output voltage and
current of the generator.
""")");

const ScriptTemplate * ScriptTemplate::Generator() {return &GeneratorScriptTemplate;}

constexpr ScriptTemplate GravitationalAccelerationScriptTemplate("GravitationalAcceleration.py", "\x01" R"(GravitationalAcceleration = """GravitationalAcceleration

g = (4 * \<pi> ** 2 * L) /
T ** 2, for small angles
""")");

const ScriptTemplate * ScriptTemplate::GravitationalAcceleration() {return &GravitationalAccelerationScriptTemplate;}

constexpr ScriptTemplate HallEffectScriptTemplate("HallEffect.py", "\x01" R"(HallEffect = """HallEffect
The Hall effect is a
phenomenon that occurs when
a current-carrying conductor
is placed in a magnetic
field. The conductor
experiences a force
perpendicular to both the
current and the magnetic
field, leading to the
creation of a potential
difference, also known as
the Hall voltage, across the
conductor.

Fe = Fb
e * E = e * v * B
E = v * B

v = I / (n * A * e)
v = I / (n * d * w * e)

V / d = B * I / (n * d * w *
e)
V = I * B / (n * w * e)

Find out if a semi-conductor
has positive or negative
charge carriers:
positive charge c
""")");

const ScriptTemplate * ScriptTemplate::HallEffect() {return &HallEffectScriptTemplate;}

constexpr ScriptTemplate HygensPrincipleScriptTemplate("HygensPrinciple.py", "\x01" R"(HygensPrinciple = """HygensPrinciple
aka propagation of waves

Each point of a wavefront is
a point source of wavelets.
These wavelets superpose and
interface to form future
wavefronts.

Explaining diffraction by
Hygen's principle:
All points on the wavefront
are sources of wavelets.
Wavelets near the edge of
the obstacle spread into the
region of geometric shadow,
creating a series minimo and
maxima as they superpose.
""")");

const ScriptTemplate * ScriptTemplate::HygensPrinciple() {return &HygensPrincipleScriptTemplate;}

constexpr ScriptTemplate InterferenceWavesScriptTemplate("InterferenceWaves.py", "\x01" R"(InterferenceWaves = """InterferenceWaves
Constructive interference
occurs when two waves are in
phase with each other and
move in the same direction.
In this scenario, the peaks
of one wave align with the
peaks of the other wave,
causing the overall
amplitude of the resulting
wave to increase. This is
because the positive parts
of the two waves add
together to create a larger
positive amplitude. For
instance, if two waves with
peak values of 1 unit each
encounter each other, the
resulting wave will have a
peak value of 2 units.
Constructive interference
occurs when the phase
difference between the waves
is an even multiple of
\<pi>(180°).

Two points in phase:
amplitudes will be added

\<delta>x - path difference
n - number of max (centre at
n = 0 ), n is natural number
\<lambda> - wavelength

\<delta>x = n * \<lambda>

Destructive interference
happens when two waves are
out of phase with each
other. One wave has a peak
where the other wave has a
trough, leading to a
decrease in the overall
amplitude of the resulting
wave. In other words, the
negative parts of the two
waves subtract from each
other, creating a smaller
overall amplitude.
Destructive interference
occurs when the phase
difference between the waves
is an odd multiple of \<pi>

Two points that are out of
phase: resultant amplitude
is the difference between
the magnitudes of the two.

n - number of minima (n = 0
for the first minima), n is
natural number

\<delta>x = (n + 0.5) *
\<lambda>

Stable interference patterns
Sources must be of the same
type of wave
Sources must be coheret
(same wavelength, freq,
phase)
Comparable amplitudes
Transverse waves must
vibrate in the same
direction (polarisation)

Two sources S1 and S2 emit
at the same frequency and
have the same phase.
Minimum: Destructive
interference
Maximum: Constructive
interference

Path difference (\<delta>x):
Difference in distance from
each source to a particular
point.

Phase difference
(\<delta>\<F ru>):
Difference in phase of the
two waves at that point.

D1(x, t) = A * sin(k * x1 -
w * t) = A * sin(\<F ru>1)
D2(x, t) = A * sin(k * x2 -
w * t) = A * sin(\<F ru>2)
\<delta>\<F ru> = \<F ru>2 -
\<F ru>1 = 
= (k * x2 - w * t) - (k * x1
- w * t)

k = 2 * \<pi> / \<lambda>
\<delta>\<F ru> = k (x2 -
x1) = k * \<delta>x
\<delta>\<F ru> = 2 * \<pi>
* \<delta>x / \<lambda>

""")");

const ScriptTemplate * ScriptTemplate::InterferenceWaves() {return &InterferenceWavesScriptTemplate;}

constexpr ScriptTemplate LenzLawScriptTemplate("LenzLaw.py", "\x01" R"(LenzLaw = """LenzLaw
Lenz’s law, in
electromagnetism, statement
that an induced electric
current flows in a direction
such that the current
opposes the change that
induced it. So this law
explains minus in the
formula: emf = - N *
\<delta>\<F ru> / \<delta>t

""")");

const ScriptTemplate * ScriptTemplate::LenzLaw() {return &LenzLawScriptTemplate;}

constexpr ScriptTemplate LightWavesScriptTemplate("LightWaves.py", "\x01" R"(LightWaves = """LightWaves
White light splits into a
range of colours called a
spectrum. The spreading
effect is known as
dispersion. Violet light has
the lowest speed in glass,
so it is refracted most.
""")");

const ScriptTemplate * ScriptTemplate::LightWaves() {return &LightWavesScriptTemplate;}

constexpr ScriptTemplate LorentzForceScriptTemplate("LorentzForce.py", "\x01" R"(LorentzForce = """LorentzForce

F = B * q * v * sin
(\<theta>)
F = n * A * l * B * B * e *
I / (n * A * e)
N = n * A * l
v = I / (n * A * e) 

Left hand rule
B - magnetic field - in the
inside side of the left hand
I - electric current (ONLY
POSITIIVE, if electron:
reversed) - in direction
with all fingers
F - Lorentz force - in
direction with big finger




""")");

const ScriptTemplate * ScriptTemplate::LorentzForce() {return &LorentzForceScriptTemplate;}

constexpr ScriptTemplate LrOrRlCircuitsScriptTemplate("LrOrRlCircuits.py", "\x01" R"(LrOrRlCircuits = """LrOrRlCircuits
A resistor–inductor circuit
(RL circuit), or RL filter
or RL network, is an
electric circuit composed of
resistors and inductors
driven by a voltage or
current source.

The current through inductor
can’t be turned on or off
instantaneously. Reason: The
change in current changes
flux, inducing an emf
opposing the change (Lenz’s
law).

the current in such a
circuit rises smoothly and
falls smoothly

We know from Lenz’s law that
inductances oppose changes
in current. There is an
alternative way to look at
this opposition that is
based on energy. Energy is
stored in a magnetic field.
It takes time to build up
energy, and it also takes
time to deplete energy;
hence, there is an
opposition to rapid change.
In an inductor, the magnetic
field is directly
proportional to current and
to the inductance of the
device. 

W or E = (L * I ** 2) / 2

Resonance frequency F0 = 1 /
sqrt(2 * pi * sqrt(L * C))
The frequency in which LC
circuit (without applied
voltage) will operate.
""")");

const ScriptTemplate * ScriptTemplate::LrOrRlCircuits() {return &LrOrRlCircuitsScriptTemplate;}

constexpr ScriptTemplate MagneticFieldsScriptTemplate("MagneticFields.py", "\x01" R"(MagneticFields = """MagneticFields
N - North - red
S - South - blue
lines go from N to S

x - line away
. - line towards (like an
arrow)


Fleming's left hand rule -
prav. lev. ruki: 
- bolshoi palez, verh -
Force experienced by the
wire
- Ukazatel., vpered - the
magnetic field
- Sred., pravo - current in
the wire 

Right hand grip rule

Radius of the curvature:
F = q * v * B
m * v ** 2 / r = q * V * B
r = m * v / (q * B)
T = 2 * \<pi> * m / (B * q)
""")");

const ScriptTemplate * ScriptTemplate::MagneticFields() {return &MagneticFieldsScriptTemplate;}

constexpr ScriptTemplate MagneticFluxScriptTemplate("MagneticFlux.py", "\x01" R"(MagneticFlux = """MagneticFlux
Magnetic flux is a
measurement of the total
magnetic field which passes
through a given area.
""")");

const ScriptTemplate * ScriptTemplate::MagneticFlux() {return &MagneticFluxScriptTemplate;}

constexpr ScriptTemplate MassSpectrometryScriptTemplate("MassSpectrometry.py", "\x01" R"(MassSpectrometry = """MassSpectrometry
There are four key stages in
the process for Mass
Spectrometry.
1. **Acceleration**: The
ions are first accelerated
until they all have the same
kinetic energy. This ensures
that all the ions entering
the velocity selector have
the same velocity.
2. **Velocity Selector**:
The ions then enter a
velocity selector. This is a
region with a uniform
electric and magnetic field.
The electric field comes
from a set of parallel
plates and is perpendicular
to the velocity of the ions,
exerting a force of `qE` on
the ions. The magnetic field
is perpendicular to both the
ion velocity and the
electric field, and its
force has a magnitude
of `qvB`, where `q` is the
charge of the ion, `v` is
the velocity, and `B` is the
strength of the magnetic
field. sin(\<theta>) = 1 ->
F _electric field_ = F
_magnetic field_ -> q * U /
d = q * v * B -> v = U / (d
* B)
3. **Deflection**: The
forces from the electric and
magnetic fields balance when
the speed of the ions
matches a certain value. For
ions with a larger speed,
the magnetic force exceeds
the electric force and those
ions are deflected out of
the beam. Similarly, ions
with a smaller speed are
deflected out of the beam in
the opposite direction,
because for those ions, the
electric force is larger
than the magnetic force.
_Uniform circular motion:_
sin(\<theta>) = 1 -> F
_centripetal_ = F _Lorentz_
-> a * m = B * q * v -> v **
2 * m / r = q * v * B -> r =
m * v / (q * B)
4. **Selection**: Therefore,
the velocity selector lets
through undeflected only
those particles with the
selected velocity. This
allows the mass spectrometer
to analyze only charged
particles with a specific
velocity.
""")");

const ScriptTemplate * ScriptTemplate::MassSpectrometry() {return &MassSpectrometryScriptTemplate;}

constexpr ScriptTemplate MotionalEmfScriptTemplate("MotionalEmf.py", "\x01" R"(MotionalEmf = """MotionalEmf
emf = - B * \<delta>A /
\<delta>t =
= - B * l * \<delta>x /
\<delta>t
v = \<delta>x / \<delta>t
emf = - B * l * v
""")");

const ScriptTemplate * ScriptTemplate::MotionalEmf() {return &MotionalEmfScriptTemplate;}

constexpr ScriptTemplate NumberOfMaximaScriptTemplate("NumberOfMaxima.py", "\x01" R"(NumberOfMaxima = """NumberOfMaxima
Theoretical limit to the
number of maxima
nmax - Highest theoretical
order of the maxima
d - distance between slits
\<lambda> - wavelength

sin(\<alpha>) <= 1
n * \<lambda> / d <= 1
nmax * \<lambda> / d = 1
nmax = d / \<lambda>

Physical limit of the
maxima:
(length of the screen and
distance to the screen)

tan(\<alpha>) = 0.5 * w / L
Determine max angle. Then
use diffr. formula
sin(\<alpha>) = n *
\<lambda> / d
nmax = sin(\<alpha>max) * d
/ \<lambda>

n - order of the maxima
Number of maxima = 2 * n + 1

DON'T FORGET ABOUTH HALFS
WHEN CALCULATING
""")");

const ScriptTemplate * ScriptTemplate::NumberOfMaxima() {return &NumberOfMaximaScriptTemplate;}

constexpr ScriptTemplate OerstedsExperimentScriptTemplate("OerstedsExperiment.py", "\x01" R"(OerstedsExperiment = """OerstedsExperiment
In 1820, Hans Christian
Oersted performed an
important experiment which
showed that there was a
connection between
electricity and
magnetism. When a current
was switched on through a
wire, it made a compass
needle turn so that it was
at right angles to the wire.

This experiment demonstrated
that a current-carrying wire
produces a magnetic
field. Oersted’s experiment
is the first which describes
this phenomenon.

Keypoints:
- An electric current
creates a circular magnetic
field as it flows through a
wire.
- The magnetic field
reverses direction when the
direction of the current is
reversed.
- The magnitude of the
current is directly
proportional to the strength
of the field.
- The strength of the
magnetic field is inversely
proportional to the distance
from the current.
- According to the
right-hand rule, the thumb
shows the direction of the
current, and four fingers
curl around the wire and
show the direction of the
magnetic field.
""")");

const ScriptTemplate * ScriptTemplate::OerstedsExperiment() {return &OerstedsExperimentScriptTemplate;}

constexpr ScriptTemplate RefractionOfLightScriptTemplate("RefractionOfLight.py", "\x01" R"(RefractionOfLight = """RefractionOfLight
c / c1 = sin i / sin i1 = n1
n1 - refractive index
""")");

const ScriptTemplate * ScriptTemplate::RefractionOfLight() {return &RefractionOfLightScriptTemplate;}

constexpr ScriptTemplate RlcCircuitScriptTemplate("RlcCircuit.py", "\x01" R"(RlcCircuit = """RlcCircuit
An RLC circuit is an
electrical circuit
consisting of a resistor
(R), an inductor (L), and a
capacitor (C), connected in
series or in parallel.

Inductor:  W = (L * I ** 2)
/ 2 Energy stored in the
magnetic field

Capacitor: W = (C * U ** 2)
/ 2 Energy stored in the
electric field 

The capacitor is charged
initially; the voltage of
this charged capacitor
causes a current to flow in
the inductor to discharge
the capacitor. Once the
capacitor is discharged, the
inductor resists any change
in the current flow, causing
the capacitor to be charged
again with the opposite
polarity. The voltage in the
capacitor eventually causes
the current flow to stop and
then flow in the opposite
direction. The result is an
oscillation, or resonance.

An important property of
this circuit is its ability
to resonate at a specific
frequency, the resonance
frequency.

Resonance occurs because
energy for this situation is
stored in two different
ways: in an electric field
as the capacitor is charged
and in a magnetic field as
current flows through the
inductor. Energy can be
transferred from one to the
other within the circuit and
this can be oscillatory. A
mechanical analogy is a
weight suspended on a spring
which will oscillate up and
down when released.

Resonance frequency: fo = 1
/ (2 * pi * sqrt(L * C))
Natural resonance angular
frequency: wo = q / sqrt(L *
C)

The frequency in which LC
circuit (no applied voltage)
will operate.
""")");

const ScriptTemplate * ScriptTemplate::RlcCircuit() {return &RlcCircuitScriptTemplate;}

constexpr ScriptTemplate SelfInductionScriptTemplate("SelfInduction.py", "\x01" R"(SelfInduction = """SelfInduction
emf = - (L * \<delta>I) / t

B = \<formula from the
booklet>

magnetic flux through one
turn: FF[Wb] = B * A

L = (N * FF) / I -> L =
\<miu>[const] * N ** 2 * A /
l[not I]

B = \<miu> * n * I -
Magnetic field in a solenoid
\<F ru> = B * A *
cos(\<theta>) - Magnetic
flux, cos(\<theta>) = 1
\<F ru> = N * \<miu> * n * I
* A = 
= \<miu> * n * \<big i> * A
* N
emf = - \<delta>\<F ru> /
\<delta>t = 
= - (\<miu> * n * A * N) *
\<big i> / t
L = - (\<miu> * A * N ** 2)
/ l

from presentation: L =
miu[const] * n ** 2 * V

In a circuit with a
solenoid: current takes
longer to reach its max.
value, light bulb brightness
only increases over time.

The larger the
self-inductance (L) of a
device, the greater its
opposition to any change.
For example: A large coil
with many turns and an iron
core has a large L, and will
not allow current to change
quickly. To avoid this, we
can counterwind coils.

The heating coils of an
electric clothes dryer can
be counter-wound so that
their magnetic fields cancel
one another, greatly
reducing the mutual
inductance with the case of
the dryer.
""")");

const ScriptTemplate * ScriptTemplate::SelfInduction() {return &SelfInductionScriptTemplate;}

constexpr ScriptTemplate SlitsWavesScriptTemplate("SlitsWaves.py", "\x01" R"(SlitsWaves = """SlitsWaves
---Single slit

a - slit width
\<lambda> - wavelength
\<alpha> - angle from centre
to first minimum

sin(\<alpha>) = \<lambda> /
a

---Two slit (Young's Double
Slit)

\<delta>x - the distance
between centre max and n-th
maximum
L - Distance between slits
and a screen
d - Distance between slits
n - order of maxima

\<delta>x = (n * L *
\<lambda>)  / d

---Multiple slits
(Diffraction grating)

d - Distance between slits
\<alpha> - the angle between
central max and n-th maximum
n - the number of maximum (n
= 0 centre) (order of
maxima)

d * sin(\<alpha>) = n *
\<lambda>

---Distance from centre to
n-th max

\<delta>y = distance between
the nth maximum and the
central max

tan(\<alpha>) = \<delta>y /
L
""")");

const ScriptTemplate * ScriptTemplate::SlitsWaves() {return &SlitsWavesScriptTemplate;}

constexpr ScriptTemplate SoundWavesScriptTemplate("SoundWaves.py", "\x01" R"(SoundWaves = """SoundWaves
Sound waves are
longitudinal.

"Loundness" the greater the
amplitude of the sound wave,
the louder the sound.
"Wavelength" is from 15 mm
to 15 m, so sound waves will
diffract round everyday
objects.
"Pitch"(Shag). The higher
the frequency, the higher
the pitch. 
""")");

const ScriptTemplate * ScriptTemplate::SoundWaves() {return &SoundWavesScriptTemplate;}

constexpr ScriptTemplate StandingWavesScriptTemplate("StandingWaves.py", "\x01" R"(StandingWaves = """StandingWaves
Form only on a string of
length L fixed at both ends
or in a pipe with both ends
closed or open if an integer
multiple of \<lambda> / 2
fits into L. (string with
two fixed nodes)

vibration is distance
between the lowest and the
highest points.
Anti-nodes is the lowest and
the highest points (peaks)

----------------------
Definition: A standing wave
is a wave that oscillates in
time but whose peak
amplitude profile does not
move in space. The peak
amplitude of the wave
oscillations at any point in
space is constant with
respect to time, and the
oscillations at different
points throughout the wave
are in phase. The locations
at which the absolute value
of the amplitude is minimum
are called nodes, and the
locations where the absolute
value of the amplitude is
maximum are called
antinodes.
Standing waves can occur due
to interference between two
waves. For instance, if a
vibrating rope tied at one
end produces a wave train
that is reflected back and
superimposed on itself as
another train of waves, the
resultant amplitude of the
two waves will be the sum of
their individual amplitudes.
At all times, there are
positions along the rope,
called nodes, at which there
is no movement at all; there
the two wave trains are
always in opposition. On
either side of a node is a
vibrating antinode. The
antinodes alternate in the
direction of displacement so
that the rope at any instant
resembles a graph of the
mathematical function called
the sine.
Standing waves can also
occur in various natural
phenomena. For example, in
the ocean, standing waves
are formed by waves with the
same wave period moving in
opposite directions. These
may form near storm centers,
or from reflection of a
swell at the shore, and are
the source of microbaroms
and microseisms.

Nodes - Positions on a
standing wave where the wave
stays in a fixed position
over time because of
destructive interference.

Atninodes - Positions on a
standing wave where the wave
vibrates with maximum
amplitude.

Internodal distance: The
internodal distance refers
to the distance between two
nodes in a structure such as
a wave. In the context of a
nerve fiber, an internodal
segment is the portion of a
nerve fiber between two
nodes of Ranvier. The term
"node" in this context
refers to a location where
the wave's amplitude is
zero, i.e., the wave has no
effect.

Internodal distance is
always \<lambda> / 2 for
standing waves. The
requirement for a standing
wave to form in these
systems is that the length
of the medium must be such
that an integer multiple of
half the wavelength fits
into that length. This
condition ensures that the
reflected waves interfere
constructively, leading to
reinforcement at certain
points (antinodes) and
destructive interference at
others (nodes), resulting in
the characteristic standing
wave pattern.

For waves where one end is
closed: - To create a
standing wave with nodes and
antinodes in this
configuration, the length of
the pipe (L) must be such
that an odd multiple of
one-fourth of the wavelength
fits into it.
Mathematically, this is
expressed as L = (2 * n + 1)
* \<lambda> / 4 , where n is
an integer. The requirement
for an odd multiple ensures
that the standing wave
pattern has an antinode at
the open end and a node at
the closed end, as the
closed end reflects the wave
with a phase change of 180
degrees (half a wavelength),
resulting in constructive
interference for odd
multiples of \<lambda>/4.

------------------

Fundamental wave is when n =
1

Standing waves with two open
ends:
L = n * \<lambda> / 2 ,
n=1,2,3...

Standing waves with one
closed and one open end:
L = n * \<lambda> / 4 ,
n=1,3,5...

How to draw:
A - antinode (i.e. waves are
far away from each other)
N - node (i.e. intersection
of waves )

two ends are open:
n = 1: ANA
n = 2: ANANA
n = 3: ANANANA

one end is closed (here left
is closed, if right then
reverse):
n = 1: NA
n = 2: NANA
n = 3: NANANA

two ends are closed:
n = 1: NAN
n = 2: NANAN
n = 3: NANANAN
""")");

const ScriptTemplate * ScriptTemplate::StandingWaves() {return &StandingWavesScriptTemplate;}

constexpr ScriptTemplate SuperpositionScriptTemplate("Superposition.py", "\x01" R"(Superposition = """Superposition
The resultant disturbance at
a point where similar waves
from several different
sources cross, is the vector
sum of individual
disturbances.
""")");

const ScriptTemplate * ScriptTemplate::Superposition() {return &SuperpositionScriptTemplate;}

constexpr ScriptTemplate TorqueScriptTemplate("Torque.py", "\x01" R"(Torque = """Torque
moment of force

\<theta> - angle between the
force and the working arm

T = F * r * sin(\<theta>)
""")");

const ScriptTemplate * ScriptTemplate::Torque() {return &TorqueScriptTemplate;}

constexpr ScriptTemplate TypesOfWavesScriptTemplate("TypesOfWaves.py", "\x01" R"(TypesOfWaves = """TypesOfWaves
Where ever there is wave
motion, there must be:
- a source of oscillation
- a material or field which
can transmit oscillations

The moving waves are called
progressive waves. Two main
types:
- Transverse waves. The
oscillations are at
right-angles to the
direction of travel. (Moving
in y axes)
- Longitudinal waves.  The
oscillations are in line
with the direction of
travel, so that a
compression ('squash') is
followed by a refraction
('stretch').
""")");

const ScriptTemplate * ScriptTemplate::TypesOfWaves() {return &TypesOfWavesScriptTemplate;}

constexpr ScriptTemplate WaveEquationScriptTemplate("WaveEquation.py", "\x01" R"(WaveEquation = """WaveEquation
y(x, t) = A * cos(2 * pi /
lambda * x {+ or -} 2 * pi /
T * t + phi)

formula like in the booklet,
but with phi (like f in ru)

A - amplitude
lambda - wavelen
T - period 
t - time
phi - shift, if functions
start neither in cos or sin

""")");

const ScriptTemplate * ScriptTemplate::WaveEquation() {return &WaveEquationScriptTemplate;}

constexpr ScriptTemplate WavesFeaturesScriptTemplate("WavesFeatures.py", "\x01" R"(WavesFeatures = """WavesFeatures
"Amplitude". This is the
magnitude (size) size of the
oscillation. (len from
x-axes to the peak).
"Frequency". This is the
number of waves emitted per
second.
"Wavelength". Distance
between one wave crest and
the next (two peaks).
"Speed". {Freq.} *
{wavelen.}

Properties of waves:
Reflection: \<aplha>1 =
\<aplha>2
Refraction: sin(\<aplha>1) /
sin(\<aplha>2) = n2 / n1 =
c1 / c2

""")");

const ScriptTemplate * ScriptTemplate::WavesFeatures() {return &WavesFeaturesScriptTemplate;}

constexpr ScriptTemplate WavesInRippleTankScriptTemplate("WavesInRippleTank.py", "\x01" R"(WavesInRippleTank = """WavesInRippleTank
Wave effects can be
investigated using a ripple
tank in which ripples travel
across the surface of
shallow water.
Reflection waves striking an
obstacle are reflected. The
angle of incidence is equal
to the angle of reflection. 
\   /
 \ / 
 ---------
left angle - angle of
incidence 
right angle - angle of
reflection
critical angle - the angle
of incidence, for which the
angle of reflection is 90.
If all light reflected -
total internal reflection.

"Refraction". When waves are
slowed down, they are
refracted (bent), provided
the angle of incidence is
not zero. In a ripple tank,
the waves can be slowed by
using a flat piece of
plastic to make the water
shallower.

"Diffraction". Waves bend
round the edges of a narrow
gap. This is called
diffraction. It is
significant if the gap size
is about a wavelength. Wider
gaps cause less diffraction.

"Interference". If two
identical sets of waves
overlap, they may either
reinforce or cancel each
other , depending on whether
they are in phase ('in
step') or out of phase.
""")");

const ScriptTemplate * ScriptTemplate::WavesInRippleTank() {return &WavesInRippleTankScriptTemplate;}

constexpr ScriptTemplate WavesGratingScriptTemplate("WavesGrating.py", "\x01" R"(WavesGrating = """WavesGrating
n * \<lambda> = d *
sin(\<alpha>)
""")");

const ScriptTemplate * ScriptTemplate::WavesGrating() {return &WavesGratingScriptTemplate;}

constexpr ScriptTemplate YoungsdoubleslitScriptTemplate("Youngsdoubleslit.py", "\x01" R"(Youngsdoubleslit = """Youngsdoubleslit
Young's double-slit
experiment is a classic
experiment in physics that
demonstrates the wave-like
behavior of light. It was
first performed by Thomas
Young in 1801, and it
involves shining light
through two closely spaced
slits onto a screen. The
resulting pattern on the
screen shows bright and dark
bands, known as interference
fringes, which reveal the
wave-like properties of
light.
The experiment begins with
sunlight being passed
through a single slit to
produce coherent light. This
light is then projected onto
a screen with two slits,
causing the light to spread
out and form interference
fringes. The pattern of
these fringes depends on the
distance between the slits,
with closer slits producing
more pronounced interference
patterns.
The interference pattern
occurs because light waves
from each slit meet at
certain points on the
screen. When the waves from
both slits are in phase
(meaning they are aligned),
they add together to create
a bright band. However, when
the waves are out of phase
(meaning one wave is leading
the other), they cancel each
other out, creating a dark
band.
This experiment provides
strong evidence for the
wave-like nature of light.

Grating and Young’s Double
Slit derivation:
Assumptions:
d << D, therefore waves from
segment from the first slit
to the maxima is parallel to
the segment from the second
slit to the maxim. Waves are
coherent because it is the
same source.

Draw a diagram where will be
two slits, \<delta>y will be
distance on the screen from
the center to the first
maxima, d is slit
separation, D is distance
between screen and slits.
\<theta> is nearest to the
slits angle of a triangle
with sides D and \<delta>y.
Small right triangle, sides:
d, n * \<lambda> (on the
L1). 

Angle between line from the
center of the slits to the
center of the screen and the
line from center of the
slits to the first maxima is
\<theta>. 
The highest angle of the
small triangle will be also
\<theta>, because angle
which is more left to
\<theta> is 90 - \<theta>,
so the highest angle of the
triangle is 90 - (90 -
\<theta>) = \<theta>.
Therefore sin(\<theta>) = n
* \<lambda> / d -> d *
sin(\<theta>) = n *
\<lambda>
tan(\<theta>) = \<delta>y /
D
we assume that angle
\<theta> is very small, so
sin(\<theta>) ==
tan(\<theta>), == is almost
equal
\<delta>y / D = n *
\<lambda> / d

Conservation of energy
Amplitude at the maximum is
double the amplitude of the
source, so intensity (I is
proportional to A ** 2) is
increased 4 times. This
increase in energy at
maximum balances the lack of
energy at minimum, so that
total energy is conserved.
""")");

const ScriptTemplate * ScriptTemplate::Youngsdoubleslit() {return &YoungsdoubleslitScriptTemplate;}

// --------------------------------------------------------------------------------------------
constexpr ScriptTemplate emptyScriptTemplate(".py", "\x01" R"(from math import *
)");

const ScriptTemplate * ScriptTemplate::Empty() {
  return &emptyScriptTemplate;
}

}
