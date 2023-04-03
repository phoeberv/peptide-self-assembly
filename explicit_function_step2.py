"""
LATEST VERSION 15 JANUARY 2023


"""
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import mdtraj as md 
#from amino_to_pdb_twochains import *
import sys

def explicit_function(PDB_input_file,final_pdb_file, DCD_output_file,temp):
    
    print ("Uploading Model...")
    pdb = PDBFile(PDB_input_file)
    #forcefield = ForceField('amber99sb.xml', 'tip3p.xml') #TIP3P model
    forcefield = ForceField('charmm36.xml', 'charmm36/water.xml')
    print ("Building Model...")
    modeller = Modeller(pdb.topology, pdb.positions) 
    
    print('Adding hydrogens...')
    modeller.addHydrogens(forcefield, pH=4.0) #Adding missing hydrogen ends
    
    print ("Adding Water Molecules...")
    modeller.addSolvent(forcefield, model = 'tip3p', padding = None) #Adding water molecules surrounding the protein strain
    
    print ("Creating System...")
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=CutoffPeriodic, nonbondedCutoff=1.0*nanometer, constraints=HBonds)
    (a,b,c) = pdb.topology.getPeriodicBoxVectors()
    system.setDefaultPeriodicBoxVectors(a, b, c)
    integrator = LangevinIntegrator(temp*kelvin, 1/picosecond, 0.001*picoseconds)  # step size: 0.001*picoseconds
    
    print("Creating Simulation Context...")
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    
    # Minimizing System 
    print("Minimizing Energy...")
    simulation.minimizeEnergy(maxIterations=25) 
    
    # Equilibrate the System -- change this
    print("Equilibrating...")
    # simulation.step(10000000) # equilibrate for 10 ns (use on the cluster)  
    simulation.step(10000000)  #equilibrate for 10 ns
    
    # Adding Reporters -- need to change this
    print("Collecting Data...")
    simulation.reporters.append(DCDReporter(DCD_output_file, 10000)) #Every 100ps record one traj data
    simulation.reporters.append(StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True, separator=' ')) 
    simulation.reporters.append(PDBReporter(final_pdb_file, 10000, enforcePeriodicBox=True)) 
    
    # Running simulation 
    simulation.step(500000000) #running for 500ns
    
def gyration(dcd_file,pdb_file):
    print("Calculating gyration radius...")
    t = md.load_dcd(dcd_file,top = pdb_file)
    gy = md.compute_rg(t)
    print (gy) #print gyration radius of protein strain
    print (t) 


if __name__ == '__main__': #explicit_function(PDB_input_file,final_pdb_file,DCD_output_file,temp)

    if len(sys.argv) == 2:
        JobName = sys.argv[1]
        (geo_list1, geo_list2) = JobName.split("_")
          
    else: 
        geo_list1 = "KGAVIL"  # change to the amino acid sequence of interest
        geo_list2 = "LIVAGK" # change the aminio acid sequence of interest
   
    
   
    #start_pdb = geo_list1 + "_" + geo_list2 + "_start.pdb"
    #amino_to_pdb_twoChains(geo_list1, geo_list2, start_pdb)
    
    
    pairsPeptides = geo_list1 + "_" + geo_list2
    start_pdb = pairsPeptides+"_start.pdb"
    final_pdb = pairsPeptides+"_final.pdb"
    traj_dcd  = pairsPeptides+"_traj.dcd"
    explicit_function(start_pdb,final_pdb,traj_dcd, temp=300)
    print("simulation complete") 
    
    #gyration('Explicit_final_25_100ns.dcd','Explicit_final_20_500.pdb') #gyration(DCD_input_file,PBD_input_file)


