#include "script_template.h"

namespace Code {

constexpr ScriptTemplate record0("B field of coil of wire (solenoid).py", "\x01" R"(aka Solenoid 
Magnitnoe pole prohodit ot + do - vnutri 
"Coil of wire", po bokam magnitnoe pole 
idet obratno i ogibaet "Coil of wire". 
)");
const ScriptTemplate * ScriptTemplate::record0_(){ return &record0;}
constexpr ScriptTemplate record1("Capacitors.py", "\x01" R"(A capacitor is a passive electronic 
component that stores electrical energy in 
an electric field by accumulating electric 
charges on two close surfaces insulated 
from each other. These two surfaces are 
often in the form of metallic plates or 
surfaces separated by a dielectric medium. 

The dielectric medium can be made of 
various materials such as glass, ceramic, 
plastic film, paper, mica, air, and oxide 
layers. The nonconducting dielectric acts 
to increase the capacitor’s charge 
capacity. 

Capacitors discharge and charge with an 
exponential decay curve. 
Charge: rost po eksponente, asimptota v 
maks zaryade 
Discharge: snizhenie po eksponente, 
asimptota v minimume 
graph voltage (y) to charge (x) shows 
stored energy under the graph 

In electrolytic capacitors, the dielectric 
is formed by the chemical action of a 
current. This gives a very thin dielectric, 
and a very high capacitance. But the 
capacitor must always be used with the same 
plate positive, or the chemical action is 
reversed. 

Capacitors have a maximum working voltage, 
above which the dielectric breaks down and 
starts to conduct. 
)");
const ScriptTemplate * ScriptTemplate::record1_(){ return &record1;}
constexpr ScriptTemplate record2("Cavendish experiment.py", "\x01" R"(The Cavendish experiment was performed in 
1797-1798 by English scientist **Henry 
Cavendish**. It was the first experiment to 
measure the force of gravity between masses 
in the laboratory, and the first to yield 
accurate values for the gravitational 
constant. 

The Cavendish experiment uses a torsion 
balance to measure the weak gravitational 
force between lead balls. A torsion balance 
consists of a bar suspended at its middle 
by a thin wire or fiber. Twisting the wire 
requires a torque that is a function of the 
wire width and material. 
The way it works is that the gravitational 
force attracting the balls together turns 
the bar, overcoming twisting resistance—or 
torque resistance—from the wire. That 
resistance is a function of angle turned 
and the torsion coefficient of the wire. At 
some angle, the torque resistance equals 
the gravitational force. 
However, the inertia of the balls causes 
them to go slightly beyond the equilibrium 
point and thus create a harmonic 
oscillation around that point. The 
oscillation is also measured by the light 
reflected from the mirror. 
The rate of oscillation is then used to 
determine the spring constant or torsion 
coefficient of the wire, which is necessary 
in the final calculation of **G**. 
)");
const ScriptTemplate * ScriptTemplate::record2_(){ return &record2;}
constexpr ScriptTemplate record3("Doppler effect.py", "\x01" R"(Doppler effect is the change in frequency 
of a wave detected by any observer because 
the wave source and the observer have 
different velocities with respect to the 
medium of the wave propagation. 

f' - shifted frequency 
v - speed of sound in medium 
vs - speed of source 
vo - speed of observer 
fs - original frequency 

Source moves towards observer: f' = v / (v 
- vs) * fs 
Source moves against observer: f' = v / (v 
+ vs) * fs 
Observer moves towards source: f' = (v + 
vo) / v * fs 
Observer moves against source: f' = (v - 
vo) / v * fs 
)");
const ScriptTemplate * ScriptTemplate::record3_(){ return &record3;}
constexpr ScriptTemplate record4("Light waves.py", "\x01" R"(White light splits into a range of colours 
called a spectrum. The spreading effect is 
known as dispersion. Violet light has the 
lowest speed in glass, so it is refracted 
most. 
)");
const ScriptTemplate * ScriptTemplate::record4_(){ return &record4;}
constexpr ScriptTemplate record5("Magnetic fields.py", "\x01" R"(N - North - red 
S - South - blue 
lines go from N to S 

x - line away 
. - line towards (like an arrow) 


Fleming's left hand rule - prav. lev. ruki: 
- bolshoi palez, verh - Force experienced 
by the wire 
- Ukazatel., vpered - the magnetic field 
- Sred., pravo - current in the wire 

Right hand grip rule 
)");
const ScriptTemplate * ScriptTemplate::record5_(){ return &record5;}
constexpr ScriptTemplate record6("Oersted's Experiment.py", "\x01" R"(In 1820, Hans Christian Oersted performed 
an important experiment which showed that 
there was a connection between electricity 
and magnetism. When a current was switched 
on through a wire, it made a compass needle 
turn so that it was at right angles to the 
wire. 

This experiment demonstrated that a 
current-carrying wire produces a magnetic 
field. Oersted’s experiment is the first 
which describes this phenomenon. 

Keypoints: 
- An electric current creates a circular 
magnetic field as it flows through a wire. 
- The magnetic field reverses direction 
when the direction of the current is 
reversed. 
- The magnitude of the current is directly 
proportional to the strength of the field. 
- The strength of the magnetic field is 
inversely proportional to the distance from 
the current. 
- According to the right-hand rule, the 
thumb shows the direction of the current, 
and four fingers curl around the wire and 
show the direction of the magnetic field. 
)");
const ScriptTemplate * ScriptTemplate::record6_(){ return &record6;}
constexpr ScriptTemplate record7("Oil-drop, Millikan exp.py", "\x01" R"(The oil-drop experiment was conducted by 
Robert Millikan and Harvey Fletcher in 1909 
in Ryerson Physical Laboratory at the 
University of Chicago. The purpose of the 
experiment was to measure the charge of a 
single electron. 

The experiment consisted of two stacked 
horizontal metal plates with insulating 
material in between. The insulating 
material had four holes, three for 
introducing light and one for microscopic 
viewing. 

A uniform electric field was created 
between the plates by applying a potential 
difference. When oil drops were sprayed, 
some of them became electrically charged 
due to friction with the nozzle. 
Alternatively, an ionizing radiation source 
like an x-ray tube could be used to charge 
the drops. 

Atomized oil drops were sprayed into a 
chamber located above the plates. Oil was 
chosen instead of water because it 
evaporates more slowly, allowing the mass 
to remain relatively constant. 

1. With the electric field turned off, the 
oil drops were introduced between the 
plates and quickly reached their terminal 
velocity due to friction with the air in 
the chamber. 
2. Upon turning on the electric field, the 
charged oil drops started to rise. This 
occurred when the electric field became 
strong enough to generate an upward 
electrical force (FE) greater than the 
gravitational force (Fg). 
3. One suitable drop was selected for 
further experimentation, and the voltage 
was alternately switched off to bring down 
the other drops while keeping the selected 
drop in the middle. 
4. The selected drop was allowed to fall at 
its terminal velocity, where the net force 
on the drop becomes zero, indicating that 
the gravitational force and the electric 
field force are equal. 

As previously mentioned, the force of 
gravity is equal to the force of the 
electric field when the drop is falling at 
its terminal velocity. 
When the voltage is turned off, the drop's 
mass is determined by its falling speed. 
Due to air resistance, the more massive 
drops fall faster than the less massive 
ones, allowing their mass to be determined 
through advanced sedimentation 
calculations. 
It was stated that the voltage (V) is 
adjusted to balance the weight of the drop 
since the applied voltage produces the 
electric field. 

Millikan achieved a precise measurement of 
the charge of the electron (qe) with an 
accuracy of 1%. Within a few years, he 
improved this measurement by a factor of 
10, arriving at a value of -1.60⋅10-19 C. 
)");
const ScriptTemplate * ScriptTemplate::record7_(){ return &record7;}
constexpr ScriptTemplate record8("Refraction of light.py", "\x01" R"(c / c1 = sin i / sin i1 = n1 
n1 - refractive index 
)");
const ScriptTemplate * ScriptTemplate::record8_(){ return &record8;}
constexpr ScriptTemplate record9("Sound waves.py", "\x01" R"(Sound waves are longitudinal. 

"Loundness" the greater the amplitude of 
the sound wave, the louder the sound. 
"Wavelength" is from 15 mm to 15 m, so 
sound waves will diffract round everyday 
objects. 
"Pitch"(Shag). The higher the frequency, 
the higher the pitch. 
)");
const ScriptTemplate * ScriptTemplate::record9_(){ return &record9;}
constexpr ScriptTemplate record10("Types of waves.py", "\x01" R"(Where ever there is wave motion, there must 
be: 
- a source of oscillation 
- a material or field which can transmit 
oscillations 

The moving waves are called progressive 
waves. Two main types: 
- Transverse waves. The oscillations are at 
right-angles to the direction of travel. 
(Moving in y axes) 
- Longitudinal waves. The oscillations are 
in line with the direction of travel, so 
that a compression ('squash') is followed 
by a refraction ('stretch'). 
)");
const ScriptTemplate * ScriptTemplate::record10_(){ return &record10;}
constexpr ScriptTemplate record11("wave equation.py", "\x01" R"(y(x, t) = A * cos(2 * pi / lambda * x {+ or 
-} 2 * pi / T * t + phi) 

formula like in the booklet, but with phi 
(like f in ru) 

A - amplitude 
lambda - wavelen 
T - period 
t - time 
phi - shift, if functions start neither in 
cos or sin 

)");
const ScriptTemplate * ScriptTemplate::record11_(){ return &record11;}
constexpr ScriptTemplate record12("Waves features.py", "\x01" R"("Amplitude". This is the magnitude (size) 
size of the oscillation. (len from x-axes 
to the peak). 
"Frequency". This is the number of waves 
emitted per second. 
"Wavelength". Distance between one wave 
crest and the next (two peaks). 
"Speed". {Freq.} * {wavelen.} 
)");
const ScriptTemplate * ScriptTemplate::record12_(){ return &record12;}
constexpr ScriptTemplate record13("Waves in ripple tank.py", "\x01" R"(Wave effects can be investigated using a 
ripple tank in which ripples travel across 
the surface of shallow water. 
Reflection waves striking an obstacle are 
reflected. The angle of incidence is equal 
to the angle of reflection. 
\ / 
\ / 
--------- 
left angle - angle of incidence 
right angle - angle of reflection 
critical angle - the angle of incidence, 
for which the angle of reflection is 90. 
If all light reflected - total internal 
reflection. 

"Refraction". When waves are slowed down, 
they are refracted (bent), provided the 
angle of incidence is not zero. In a ripple 
tank, the waves can be slowed by using a 
flat piece of plastic to make the water 
shallower. 

"Diffraction". Waves bend round the edges 
of a narrow gap. This is called 
diffraction. It is significant if the gap 
size is about a wavelength. Wider gaps 
cause less diffraction. 

"Interference". If two identical sets of 
waves overlap, they may either reinforce or 
cancel each other , depending on whether 
they are in phase ('in step') or out of 
phase. 
)");
const ScriptTemplate * ScriptTemplate::record13_(){ return &record13;}



constexpr ScriptTemplate emptyScriptTemplate(".py", "\x01" R"(from math import *
)");

constexpr ScriptTemplate squaresScriptTemplate("bibbob.py", "\x01" R"(here will be a lot of
interesting
and useful


text)");

constexpr ScriptTemplate mandelbrotScriptTemplate("mandelbrot.py", "\x01" R"(# This script draws a Mandelbrot fractal set
# N_iteration: degree of precision
import kandinsky
def mandelbrot(N_iteration):
  for x in range(320):
    for y in range(222):
# Compute the mandelbrot sequence for the point c = (c_r, c_i) with start value z = (z_r, z_i)
      z = complex(0,0)
# Rescale to fit the drawing screen 320x222
      c = complex(3.5*x/319-2.5, -2.5*y/221+1.25)
      i = 0
      while (i < N_iteration) and abs(z) < 2:
        i = i + 1
        z = z*z+c
# Choose the color of the dot from the Mandelbrot sequence
      rgb = int(255*i/N_iteration)
      col = kandinsky.color(int(rgb),int(rgb*0.75),int(rgb*0.25))
# Draw a pixel colored in 'col' at position (x,y)
      kandinsky.set_pixel(x,y,col))");

constexpr ScriptTemplate polynomialScriptTemplate("polynomial.py", "\x01" R"(from math import *
# roots(a,b,c) computes the solutions of the equation a*x**2+b*x+c=0
def roots(a,b,c):
  delta = b*b-4*a*c
  if delta == 0:
    return -b/(2*a)
  elif delta > 0:
    x_1 = (-b-sqrt(delta))/(2*a)
    x_2 = (-b+sqrt(delta))/(2*a)
    return x_1, x_2
  else:
    return None)");

constexpr ScriptTemplate parabolaScriptTemplate("parabola.py", "\x01" R"(from matplotlib.pyplot import *
from math import *

g=9.81

def x(t,v_0,alpha):
  return v_0*cos(alpha)*t
def y(t,v_0,alpha,h_0):
  return -0.5*g*t**2+v_0*sin(alpha)*t+h_0

def vx(v_0,alpha):
  return v_0*cos(alpha)
def vy(t,v_0,alpha):
  return -g*t+v_0*sin(alpha)

def t_max(v_0,alpha,h_0):
  return (v_0*sin(alpha)+sqrt((v_0**2)*(sin(alpha)**2)+2*g*h_0))/g

def simulation(v_0=15,alpha=pi/4,h_0=2):
  tMax=t_max(v_0,alpha,h_0)
  accuracy=1/10**(floor(log10(tMax))-1)
  T_MAX=floor(tMax*accuracy)+1
  X=[x(t/accuracy,v_0,alpha) for t in range(T_MAX)]
  Y=[y(t/accuracy,v_0,alpha,h_0) for t in range(T_MAX)]
  VX=[vx(v_0,alpha) for t in range(T_MAX)]
  VY=[vy(t/accuracy,v_0,alpha) for t in range(T_MAX)]
  for i in range(T_MAX):
    arrow(X[i],Y[i],VX[i]/accuracy,VY[i]/accuracy)
  grid()
  show())");

const ScriptTemplate * ScriptTemplate::Empty() {
  return &emptyScriptTemplate;
}

const ScriptTemplate * ScriptTemplate::Squares() {
  return &squaresScriptTemplate;
}

const ScriptTemplate * ScriptTemplate::Mandelbrot() {
  return &mandelbrotScriptTemplate;
}

const ScriptTemplate * ScriptTemplate::Polynomial() {
  return &polynomialScriptTemplate;
}

const ScriptTemplate * ScriptTemplate::Parabola() {
  return &parabolaScriptTemplate;
}

}
