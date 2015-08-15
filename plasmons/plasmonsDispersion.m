# not a function file
1;

#{
  This program evaluates the dispersion curve for the 
  optical plasmons propagating in the topological insulator
  surface. The topological insulator is surrounded by two 
  dielectric media with dielectric constants epsilonT and 
  epsilonB.
  The calculations are based on the formula given in the paper
  T. Stauber and G. Gomez-Santos and L. Brey
  "Spin-charge separation of plasmonic excitations in 
  thin film topological insulators"
  Physical Review B 88, 205427 (2013)
  
#}


# fixed constants-------------------------------
global Hbar = 1.05457173E-34; # reduced Planck Constant
global Pi = 3.14159265358979323846; # no comment
global Epsilon0 = 8.854187817E-12; # vacuum permittivity 
global Ec = 1.60217657E-19; # electron charge magnitude
#------------------------------------------------


# Experimental parameters-------------------------
global vF = 6.0E5; # m/s, Fermi velocity
global kF = 7.5E8#7.5E8##1.37E9#7.5E8; # 1/m, Fermi wave-vector
global mu = Hbar*vF*kF;#0.296*Ec; # Chemical potential for Top and Bottom surfaces
global epsilonT = 1;#6.7#5.76; # dielectric constant for Top layer
global epsilonB = 10.0; # dielectric constant for Bottom layer
global epsilonTI = 25.0; # dielectric constant for Topological Insulator
global d = 15E-9; # Topological Insulator film thickness
#--------------------------------------------------

global alpha = 7.29 #fine structure of a Dirac system e^2/(eps0*h*vF)

function y = G( x ) 
# auxilary function, needed for further evaluation
# of polarizability 
  if x < 1
    y = -1;
  else
    y = x*sqrt(x^2-1) - 1/cosh(x);
  endif
endfunction

function y = F( q_, w_ )
# another auxilary function of dimensionless wavenumber and 
# frequency
  y  = 1 - (q_^2) /(8*sqrt(w_^2-q_^2))*( G( (2+w_)/q_ ) - G( (2-w_)/q_ ) );
endfunction

function y = N( q_ )
# auxilary function, needed for further evaluation
# of Coulomb interaction
  global epsilonB;
  global epsilonT;
  global epsilonTI;
  global d;
  global mu;
  global Hbar;
  global vF;
  d_ = d*mu/(Hbar*vF); # dimensionless thickness of the film
  y = epsilonTI*(epsilonT+epsilonB)*cosh(q_*d_) + (epsilonT*epsilonB+epsilonTI^2)*sinh(q_*d_);
endfunction


function y = Beta( q_, w_ )
  global alpha;
  global epsilonTI;
  y = -alpha*epsilonTI*F(q_, w_)/(q_*N(q_));
endfunction


function y = C( q_ )
  global epsilonB;
  global epsilonTI;
  global d;
  global mu;
  global Hbar;
  global vF;
  d_ = d*mu/(Hbar*vF); # dimensionless thickness of the film
  y = cosh(q_*d_)+epsilonB/epsilonTI*sinh(q_*d_);
endfunction

function y = D( q_ )
  global epsilonT;
  global epsilonTI;
  global d;
  global mu;
  global Hbar;
  global vF;
  d_ = d*mu/(Hbar*vF); # dimensionless thickness of the film
  y = cosh(q_*d_)+epsilonT/epsilonTI*sinh(q_*d_);
endfunction

function y = H( q_, w_ )
  y = (1 - C(q_)*Beta(q_, w_))*(1 - D(q_)*Beta(q_, w_)) - Beta(q_, w_)*Beta(q_, w_);
endfunction

function y = wZeroApprox( q_ )
  # normalized plasmonic frequency in the zero order approximation;
  # See formula (9) in the reference
  global alpha;
  global epsilonT;  
  global epsilonB;
  y = sqrt( alpha*q_/( epsilonB+epsilonT) );
endfunction

function y = wFirstApprox( q_ )
  # normalized plasmonic frequency in the first order approximation
  # See formula (12) in the reference
  global alpha;
  global epsilonT;
  global epsilonTI;
  global epsilonB;
  global d;
  global mu;
  global Hbar;
  global vF;
  d_ = d*mu/(Hbar*vF); # dimensionless thickness of the film
  y = sqrt( alpha*q_/( epsilonB+epsilonT + epsilonTI*q_*d_ ) );
endfunction

# solving for the dispersion relation
# step 1: defining the wave-vectors 
numberOfPoints = 400;
qOverkF = 0.5;
Q = linspace(1E-4*qOverkF, qOverkF, numberOfPoints);
# step 2: preparing frequency vectors for
# W: exact value, from solving the equation (7)
# WZero: very approximate value, formula (9) in the reference
# WFirst: approximate value, formula (12) in the reference
W = zeros(1, numberOfPoints);
WZero = zeros(1, numberOfPoints);
WFirst = zeros(1, numberOfPoints);
# step 3: solving for zeros of determinant for every wave-number
for k = 1:numberOfPoints
  q = Q(k);
  # evaluating approximate values for the plasmonic frequencies.
  WZero(k) = wZeroApprox(q);
  WFirst(k) = wFirstApprox(q);  
  # defining temporary anonymous function of one variable  
  f = @(x) H( q, x );  
  #[x, fval, info] = fsolve( f, [wOpt0(k)] );
  guess = (WZero(k) + WFirst(k))/2;  
  [x, fval, info, output] = fsolve( f, guess );
  W(k) = x;
endfor

plot(Q, W, "color", "red", "linewidth", 2, "linestyle", "-");
xlabel("q/kF");
ylabel("hw/mu");
hold on;
plot(Q, WFirst, "color", "green", "linewidth", 2, "linestyle", "-.");
hold on;
plot(Q, WZero, "color", "blue", "linewidth", 2, "linestyle", "--");
