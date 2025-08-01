
# ******************************************************************FOR H2 MOLECULE**********************************************************************************************************
import pennylane as qml
from pennylane import numpy as np

# Define the H2 molecule
symbols = ["H", "H"]
coordinates = np.array([0.0, 0.0, -0.6614, 0.0, 0.0, 0.6614])

charge = 0
multiplicity = 1  # H2 has a singlet ground state
basis_set = "sto-3g"  # Minimal basis set

# Create the Molecule object
molecule = qml.qchem.Molecule(
    symbols,
    coordinates,
    charge=charge,
    mult=multiplicity,
    basis_name=basis_set
)

# Define the number of active electrons and orbitals
active_electrons = 2  # 2 electrons in H2
active_orbitals = 2  # 2 orbitals, each hydrogen contributes a 1s orbital

# Compute the Hamiltonian and number of qubits
H, qubits = qml.qchem.molecular_hamiltonian(
    molecule,
    active_electrons=active_electrons,
    active_orbitals=active_orbitals,
)

print("Number of qubits required to perform quantum simulations: {:}".format(qubits))
print("Hamiltonian of the H₂ molecule:")
print(H)

# Generate the Hartree-Fock state
hf_state = qml.qchem.hf_state(active_electrons, qubits)

# Generate single and double excitations
singles, doubles = qml.qchem.excitations(active_electrons, qubits)

print("Hartree-Fock state:")
print(hf_state)
print("Single excitations:")
print(singles)
print("Double excitations:")
print(doubles)

# Map excitations to the wires for the UCCSD circuit
s_wires, d_wires = qml.qchem.excitations_to_wires(singles, doubles)

# Define the device
dev = qml.device("default.qubit", wires=qubits)

# Define the qnode and ansatz
@qml.qnode(dev)
def circuit(params, wires, s_wires, d_wires, hf_state):
    qml.UCCSD(params, wires, s_wires, d_wires, hf_state)
    return qml.expval(H)

# Define the initial values of the circuit parameters
params = np.zeros(len(singles) + len(doubles))

# Define the optimizer
optimizer = qml.GradientDescentOptimizer(stepsize=0.5)

# Optimize the circuit parameters and compute the energy
for n in range(21):  # Increase the number of iterations for more accurate energy
    params, energy = optimizer.step_and_cost(circuit, params,
    wires=range(qubits), s_wires=s_wires, d_wires=d_wires, hf_state=hf_state)
    if n % 2 == 0:
        print("step = {:},  E = {:.8f} Ha".format(n, energy))
      
#******************************************************************FOR O2 MOLECULE**********************************************************************************************************
import pennylane as qml
from pennylane import numpy as np
from pennylane import qchem
import openfermion
import openfermionpyscf

# Define the molecule
symbols = ["O", "O"]
coordinates = np.array([[0.0, 0.0, 0.0], [2.282, 0.0, 0.0]])  # Approximate bond length in bohr atomic unit

charge = 0
multiplicity = 3 # O₂ molecule has a triplet ground state
basis_set = "sto-3g"

molecule = qml.qchem.Molecule(
    symbols,
    coordinates,
    charge=charge,
    mult=multiplicity,
    basis_name=basis_set
)

# Define the number of active electrons and orbitals
active_electrons =12 
active_orbitals = 8  

# Compute the Hamiltonian and number of qubits
H, qubits = qml.qchem.molecular_hamiltonian(
    molecule,
    active_electrons=active_electrons,
    active_orbitals=active_orbitals,
    method = 'openfermion'
)

print(H)
print(qubits)

# Generate the Hartree-Fock state
hf_state = qml.qchem.hf_state(active_electrons, qubits)

# Generate single and double excitations
singles, doubles = qml.qchem.excitations(active_electrons, qubits)

print(hf_state)
print(singles)
print(doubles)


# Map excitations to the wires for the UCCSD circuit
s_wires, d_wires = qml.qchem.excitations_to_wires(singles, doubles)

# Define the device
dev = qml.device("lightning.qubit", wires=qubits)

# Define the qnode and ansatz
@qml.qnode(dev)
def circuit(params, wires, s_wires, d_wires, hf_state):
    qml.UCCSD(params, wires, s_wires, d_wires, hf_state)
    return qml.expval(H)

# Define the initial values of the circuit parameters
params = np.zeros(len(singles) + len(doubles))

# Define the optimizer
optimizer = qml.GradientDescentOptimizer(stepsize=0.5)


# Optimize the circuit parameters and compute the energy
for n in range(21):  # Increase the number of iterations for more accurate energy
    params, energy = optimizer.step_and_cost(circuit, params,
    wires=range(qubits), s_wires=s_wires, d_wires=d_wires, hf_state=hf_state)
    if n % 2 == 0:
        print("step = {:},  E = {:.8f} Ha".format(n, energy))


#******************************************************************FOR H2O MOLECULE**********************************************************************************************************
import pennylane as qml
from pennylane import numpy as np

symbols = ["H", "O", "H"]
coordinates = np.array([

[1.1103, 1.4235, 0.0], # Hydrogen atom 1
    [0.0, 0.0, 0.0], # Oxygen atom
[1.1103, -1.4235, 0.0] # Hydrogen atom 2
])

charge = 0
multiplicity = 1
basis_set = "sto-3g"

molecule = qml.qchem.Molecule(
    symbols,
    coordinates,
    charge=charge,
    mult=multiplicity,
    basis_name=basis_set
)

H, qubits = qml.qchem.molecular_hamiltonian(
    molecule,
    active_electrons=4
    active_orbitals=4
)

print("Number of qubits required to perform quantum simulations: {:}".format(qubits))
print("Hamiltonian of the water molecule")
print(H)

electrons = 4
# Define the HF state
hf_state = qml.qchem.hf_state(electrons, qubits)

# Generate single and double excitations
singles, doubles = qml.qchem.excitations(electrons, qubits)

print(hf_state)
print(singles)
print(doubles)

# Map excitations to the wires the UCCSD circuit will act on
s_wires, d_wires = qml.qchem.excitations_to_wires(singles, doubles)

# Define the device
dev = qml.device("default.qubit", wires=qubits)

# Define the qnode and ansatz
@qml.qnode(dev)
def circuit(params, wires, s_wires, d_wires, hf_state):
    qml.UCCSD(params, wires, s_wires, d_wires, hf_state)
    return qml.expval(H)

# Define the initial values of the circuit parameters
params = np.zeros(len(singles) + len(doubles))

# Define the optimizer
optimizer = qml.GradientDescentOptimizer(stepsize=0.5)

# Optimize the circuit parameters and compute the energy
for n in range(21):      #CAN INCREASE THE NUMBER OF ITERATIONS FOR MORE ACCURATE ENERGY
    params, energy = optimizer.step_and_cost(circuit, params,
    wires=range(qubits), s_wires=s_wires, d_wires=d_wires, hf_state=hf_state)
    if n % 2 == 0:
        print("step = {:},  E = {:.8f} Ha".format(n, energy))
      
