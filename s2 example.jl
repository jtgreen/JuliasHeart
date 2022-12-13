include("Tor_Ord_Land.jl")

const tstop1 = [5., 300.]
const tstop2 = [6., 301.]


# starts
function condition(u,t,integrator)
    t in tstop1
end

# stops
function condition2(u,t,integrator)
    t in tstop2
end

function affect!(integrator)
    integrator.u[50] = -53.0
end
  
function affect2!(integrator)
    integrator.u[50] = 0.0
end

save_positions = (true,true)

cb = DiscreteCallback(condition, affect!, save_positions=save_positions)

save_positions = (false,true)

cb2 = DiscreteCallback(condition2, affect2!, save_positions=save_positions)

cbs = CallbackSet(cb,cb2)


# % Run for ToR-ORd+Land electro-mechanical model
# % (Margara, F., Wang, Z.J., Levrero-Florencio, F., Santiago, A., Vázquez, M., Bueno-Orovio, A.,
# % and Rodriguez, B. (2021). In-silico human electro-mechanical ventricular modelling and simulation for
# % drug-induced pro-arrhythmia and inotropic risk assessment. Progress in Biophysics and Molecular Biology).
# % https://doi.org/10.1016/j.pbiomolbio.2020.06.007
# %% Setting parameters

param = Dict(	"bcl"   =>      1000, 					# % basic cycle length in ms
            "model"     =>      model_ToRORd_Land, 	    # % which model is to be used
            "verbose"   =>      true, 				    # % printing numbers of beats simulated.
            "cellType"  =>      0, 					    # %0 endo, 1 epi, 2 mid
            "callback"  =>      cbs, 				    # % callback for pacing
            "tstops"    =>      [5.;6.;300.;301.], 			    # % stop simulation after 10000 ms
)

options = [] # % parameters for ode15s - usually empty
beats = 5 # % number of beats JTG: (200 default in original code)
ignoreFirst = beats - 1 # % this many beats at the start of the simulations are ignored when extracting the structure of simulation outputs (i.e., beats - 1 keeps the last beat).
X0 = getStartingState("m_endo") # % starting state - can be also m_mid or m_epi for midmyocardial or epicardial cells respectively.

# % Simulation and extraction of outputs
X = modelRunner(X0, options, param, beats, ignoreFirst)
p = plot(X, vars=1)
display(p)
# currents = getCurrentsStructure(time, X, beats, param, 0)
# % ActiveTension = currents.Ta