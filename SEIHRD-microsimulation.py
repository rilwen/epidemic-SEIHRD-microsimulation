# -*- coding: utf-8 -*-
"""
SEIHR microsimulation.

"""

import numpy as np
import random as rd
import enum
import time

# Maximum hospital capacity above which the fatality rate increases because of
# overburdening of hospitals (downscaled to the simulated population size).
NHS_OVERLOAD = 2000

# List of numbers of families of a given size.
# The i-th element is the number of families of size i+1.
FAMILIES = [int(1e6),int(1.3e6),int(0.2e6)]


class Status(enum.IntEnum):
    """SEIHRD stages."""
    
    # Susceptible to infection through exposure to an infected person with
    # whom they were in contact.
    SUSCEPTIBLE = 0
    
    # Has been in contact with an infected person on a given day
    # (but isn't spreading).
    EXPOSED = 1
    
    # Developed the infection following the exposure and virus incubation period
    # (is spreading).
    INFECTED = 2
    
    # Infected and admitted to hospital.
    HOSPITALISED = 3
    
    #Recovered from infection and immune.
    RECOVERED = 4
    
    # Deceased at hospital (only hospital cases can die).
    DEAD = 5



def simulate(n_days, families, n_init_exposed,max_n_contacts,max_freq):
    """Simulates SEIHRD model.
    
    Calculates consecutive states of the simulation for each person based on
    probabilities of transitions between SEIHRD states.
    
    Args:
        n_days: Length of simulation in days.
        families: List whose i-th element is the number of families of size i+1.
        n_init_exposed: Initial number of persons exposed to the virus.
    
    Returns:
        state: Numpy array with shape (n_days, number of persons) containing
            the simulated trajectories as persons' SEIHRD integers states.
    """
    #t0=time.time()
    
    trajectory_log = open(
            "simulate_contacts%d_freq%d.csv"%(max_n_contacts,max_freq), "w")
    
    # Number of persons in the simulation.
    n_persons= count_persons(families)
        
    # Microsimulation state with simulated days in rows and persons in columns.
    state = np.zeros((n_days, n_persons), dtype=int)
    
    # Randomly expose to virus the initial number of persons (n_init_exposed).
    state[0,0:n_init_exposed] = Status.EXPOSED
    rd.shuffle(state[0,:])
    
    print("#0: Assign family ID and add famility members to contacts")
    persons_familyID, persons_contacts = assign_family_ID_and_members(families)
    
    print("#1: Add external contacts")
    persons_contacts = assign_external_contacts(
            persons_contacts,persons_familyID,max_n_contacts,max_freq)
    
    print("#2: Simulate")
    np.savetxt(trajectory_log, state[0])
    
    # SEIHRD model parameters:
    min_HR_time = 14 # minimum time at hospital before recovery
    min_HD_time = 5 # minimum time at hospital before death
    min_IR_time = 7 # minimum recovery time (in days)
    EI_time = 4 # virus incubation time (in days)
    EI_prob = 0.5 # probability of developing infection after incubation time
    ES_prob = 1-EI_prob # p. of not developing infection after incubation time
    HR = 1/35 # daily probability of recovery at hospital
    IR = 1/21 # daily probability of recovery without hospitalisation
    IH = 1/200 # daily probability of an infected person going to hospital
    base_HD = HR*0.16 # daily probability of death at hospital
    
    for day in range(1,n_days):
        print("Day ",day)
        n_hospitalised = sum(state[day-1]==Status.HOSPITALISED)
        
        if n_hospitalised < NHS_OVERLOAD:
            HD = base_HD
        else:
            HD = 3*base_HD
        
        for i in range(n_persons):
            if state[day-1,i] in (Status.DEAD,Status.RECOVERED):
                # DEAD and RECOVERED are terminal states.
                state[day,i] = state[day-1,i]
        
            elif state[day-1,i] in (Status.EXPOSED, Status.SUSCEPTIBLE):
                
                # Calculate the probability of getting exposed on this day.
                exposure_frequencies = [freq for c,freq in persons_contacts[i]
                    if state[day-1,c]==Status.INFECTED]
                exposure_probability = min(1,sum(exposure_frequencies))
                
                if state[day-EI_time,i] != Status.EXPOSED:
                    state[day,i] = rd.choices([Status.EXPOSED,Status.SUSCEPTIBLE],
                         weights=[exposure_probability,1-exposure_probability])[0]                
                else:
                    # Once the virus incubation time has passed, the person is
                    # at risk of infection. If they don't become infected, they
                    # become susceptible or can be exposed again.
                    state[day,i] = rd.choices([Status.INFECTED,Status.SUSCEPTIBLE,Status.EXPOSED],
                        weights=[EI_prob,(1-exposure_probability)*ES_prob, exposure_probability*ES_prob])[0]

            elif state[day-1,i] == Status.INFECTED:
                # If minimum infection time has passed the person can recover with immunity.
                if day >= min_IR_time and state[day-min_IR_time,i] == Status.INFECTED:
                    IR_prob = IR
                else:
                    IR_prob = 0
                IH_prob = IH
                II_prob = 1 - (IR_prob + IH_prob)
                state[day,i] = rd.choices([Status.RECOVERED,Status.HOSPITALISED,Status.INFECTED],
                     weights=[IR_prob,IH_prob,II_prob])[0]

            elif state[day-1,i] == Status.HOSPITALISED:
                
                # If minimum hospitalisation time has passed the person can recover with immunity.
                if day >= min_HR_time and state[day-min_HR_time,i] == Status.HOSPITALISED:
                    HR_prob = HR
                else:
                    HR_prob = 0
                    
                # If minimum hospitalisation time has passed the person can die.
                if day >= min_HD_time and state[day-min_HD_time,i] == Status.HOSPITALISED:
                    HD_prob = HD
                else:
                    HD_prob = 0                    
                HH_prob = 1-(HR_prob+HD_prob)
                state[day,i] = rd.choices([Status.HOSPITALISED,Status.RECOVERED,Status.DEAD],weights=[HH_prob,HR_prob,HD_prob])[0]
                
        np.savetxt(trajectory_log, state[day]) # Save this day's states.

    #tf=time.time()
    #print("Simulation time [s]", tf-t0)
    trajectory_log.close()
    return state


def count_persons(families):
    """
    Calculates the number of persons in the simulation.
    Args:
        families: List whose i-th element is the number of families of size i+1.
    Returns:
        n_persons: Total number of persons in the simulation.
    """
    n_persons = sum([(i+1)*n for i,n in enumerate(families)])
    return n_persons


def assign_family_ID_and_members(families):
    """
    Assigns the same unique index to family members and adds family members to list of contacts.
    Args:
        families: List whose i-th element is the number of families of size i+1.
    Returns:
        persons_familyID: List of persons' family IDs.
        persons_contacts: List of lists of persons' contacts and daily contact frequencies.
    """
    n_persons = count_persons(families)
    persons_familyID = np.zeros(n_persons, dtype=int)
    persons_contacts = [[] for _ in range(n_persons)]
    familyID = 0
    personID = 0
    for i in range(len(families)):
        family_size = i+1
        for n in range(families[i]):
            persons_familyID[personID:personID+family_size] = familyID
            # Assume infinite contact frequency with family members.
            family_contacts = [(p,np.inf) for p in range(personID,personID+family_size)]
            for p in range(personID,personID+family_size):
                persons_contacts[p] = [c for c in family_contacts if c[0]!=p]
            personID += family_size
            familyID +=1
    return persons_familyID, persons_contacts


def assign_external_contacts(persons_contacts, persons_familyID,max_n_contacts,max_freq):
    """
    Adds randomly selected external contacts, i.e. outside family, to each person's
    contact list. Each person initiates a random number of external contacts with
    randomly chosen persons. Each contact has a fixed frequency chosen
    randomly from a defined range.
    Args:
        persons_familyID: List of persons' family IDs.
        persons_contacts: Current list of lists of persons' contacts and daily contact frequencies.
        max_n_contacts: Maximum number of external contacts initiated by the person.
        max_freq: Maximum weekly frequency for each external contact
    Returns:
        persons_contacts: List of lists of persons' contacts and daily frequencies
        extended by external contacts. persons_contacts[i][k] is a tuple (p2, f)
        meaning that person p2 is the k-th contact of person p with frequency f.
    """
    n_persons = len(persons_contacts)    
    for p in range(n_persons):
        n_contacts = rd.randint(0,max_n_contacts) # Hence, the average number of external contacts (initiated + received) per person is max_n_contacts
        c = 0
        contacts_set = set()
        while c < n_contacts:
            p2 = rd.randint(0,n_persons-1)
            if persons_familyID[p2] != persons_familyID[p] and p2 not in contacts_set:
                contacts_set.add(p2)
                frequency = rd.randint(1,max_freq)/7
                # Initiate the contact from p to p2.
                persons_contacts[p].append((p2,frequency))
                # p2 reciprocates.
                persons_contacts[p2].append((p,frequency))
                c +=1
    return persons_contacts
        
def count_stock(state,person_status):
    """ Counts the number of person with given status on each day. """
    return np.sum(state==person_status,1)
    

simulation_nsi = simulate(365,FAMILIES,100,5,5) # no self-isolation
simulation_si1 = simulate(365,FAMILIES,100,3,5) # light self-isolation
simulation_si2 = simulate(365,FAMILIES,100,3,1) # moderate self-isolation
simulation_si3 = simulate(365,FAMILIES,100,2,1) # harsh self-isolation
simulation_si4 = simulate(1000,FAMILIES,100,1,1) # extreme self-isolation


import matplotlib.pyplot as plt
S = count_stock(simulation_nsi ,Status.SUSCEPTIBLE)
E = count_stock(simulation_nsi ,Status.EXPOSED)
I = count_stock(simulation_nsi ,Status.INFECTED)
H = count_stock(simulation_nsi ,Status.HOSPITALISED)
R = count_stock(simulation_nsi ,Status.RECOVERED)
D = count_stock(simulation_nsi ,Status.DEAD)


plt.figure(figsize = (8,6))
#plt.plot(S,label='Susceptible')
#plt.plot(E,label='Exposed')
#plt.plot(I,label='Infected')
#plt.plot(R,label='Recovered')
#plt.plot(H,label='Hospitalised')
plt.plot(D, label='Dead')
plt.legend()
plt.grid()
plt.xlabel("days")