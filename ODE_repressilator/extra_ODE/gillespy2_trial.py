# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 09:35:50 2020

@author: pablo
"""

#This is just a trial, we did not manage to get it working
#the equations in this model are from: Transient Dynamics of Genetic Regulatory Networks (Bennet 2007)

import gillespy2
import numpy 



class Repressilator(gillespy2.Model):
    def __init__(self, parameter_values=None):
        # First call the gillespy2.Model initializer.
        gillespy2.Model.__init__(self, name='Repressilator')

        # Define parameters for the rates of creation and dissociation.
        alpha=gillespy2.Parameter(name="alpha", expression=216)
        alpha_0=gillespy2.Parameter(name="alpha_0",expression=0.216)
        beta=gillespy2.Parameter(name="beta",expression=0.8)
        self.add_parameter([alpha, alpha_0, beta])

        # Define variables for the molecular species representing M and D.
        m_lacI= gillespy2.Species(name='m_lacI', initial_value=30)
        p_lacI= gillespy2.Species(name='p_lacI', initial_value=0)
        m_tetR= gillespy2.Species(name='m_tetR', initial_value=0)
        p_tetR= gillespy2.Species(name='p_tetR', initial_value=0)
        m_cI= gillespy2.Species(name='m_cI', initial_value=0)
        p_cI= gillespy2.Species(name='p_cI', initial_value=0)
        self.add_species([m_lacI, p_lacI, p_tetR, m_tetR, m_cI, p_cI])

        # The list of reactants and products for a Reaction object are each a
        # Python dictionary in which the dictionary keys are Species objects
        # and the values are stoichiometries of the species in the reaction.
        
        m_1 = gillespy2.Reaction(name="m_1", rate=alpha/(1 + p_cI), reactants={m_lacI:1}, products={m_lacI:2+alpha_0})
        m_2 = gillespy2.Reaction(name="m_2", rate=1, reactants={m_lacI:1}, products={m_lacI:0})
        p_1 = gillespy2.Reaction(name="p_1", rate=beta, reactants={p_lacI:1}, products={p_lacI:2+m_lacI})                       
        p_2 = gillespy2.Reaction(name="p_2", rate=beta, reactants={p_lacI:1}, products={p_lacI:0})                                   
        
        m_3 = gillespy2.Reaction(name="m_3", rate=alpha/(1+p_lacI), reactants={m_tetR:1}, products={m_tetR:2+alpha_0})
        m_4 = gillespy2.Reaction(name="m_4", rate=1, reactants={m_tetR:1}, products={m_tetR:0})
        p_3 = gillespy2.Reaction(name="p_3", rate=beta, reactants={p_tetR:1}, products={p_tetR:2+m_tetR})                  
        p_4 = gillespy2.Reaction(name="p_4", rate=beta, reactants={p_tetR:1}, products={p_tetR:0})
        
        m_5 = gillespy2.Reaction(name="m_5", rate=alpha/(1+p_tetR), reactants={m_cI:1}, products={m_cI:2+alpha_0})
        m_6 = gillespy2.Reaction(name="m_6", rate=1, reactants={m_cI:1}, products={m_cI:0})
        p_5 = gillespy2.Reaction(name="p_5", rate=beta, reactants={p_cI:1}, products={p_cI:2+m_cI})                  
        p_6 = gillespy2.Reaction(name="p_6", rate=beta, reactants={p_cI:1}, products={p_cI:0})
        
        self.add_reaction([m_1,m_2,m_3,m_4,m_5,m_6,p_1,p_2,p_3,p_4,p_5,p_6])
        
        # Set the timespan for the simulation.
        self.timespan(numpy.linspace(0, 100, 101))
       
model = Repressilator()
results = model.run(number_of_trajectories=10)
results.plot()