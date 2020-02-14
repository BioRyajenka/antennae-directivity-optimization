import random
from deap import base, creator, tools
import numpy as np
from math import sqrt

creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

toolbox = base.Toolbox()
#toolbox.register("individual", tools.initIterate, creator.Individual, createIndividual)
#toolbox.register("population", tools.initRepeat, list, toolbox.individual)

def update_fitnesses(population):
    invalid_ind = [ind for ind in population if not ind.fitness.valid]
    fitnesses = map(toolbox.evaluate, invalid_ind)
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = (fit,)

def get_ind_function(limits):
    def ind_function():
        return list(map(lambda t: random.uniform(t[0], t[1]), limits))
    return ind_function

def get_mutate_func(limits):
    def mutate_func(ind):
        for i in range(len(ind)):
            mutated_value = np.random.normal(ind[i], (limits[i][1] - limits[i][0]) / 6.0 / sqrt(len(ind)))

            #print("hmm {} {} ! {} {} {}".format(ind[i], mutated_value, limits[i][0], limits[i][1], (limits[i][1] - limits[i][0]) / 6.0))

            mutated_value = limits[i][1] if mutated_value > limits[i][1] else mutated_value
            mutated_value = limits[i][0] if mutated_value < limits[i][0] else mutated_value

            #print("hmm2 {}".format(mutated_value))
            assert mutated_value >= limits[i][0] and mutated_value <= limits[i][1]
            ind[i] = mutated_value

        return ind
    return mutate_func

def one_plus_one_onecall(evaluate_func, limits, max_plato_time=100, log=False, callback=None):
    toolbox.register("individual", tools.initIterate, creator.Individual, get_ind_function(limits))
    toolbox.register("evaluate", evaluate_func)
    toolbox.register("select", tools.selBest, k=1)
    toolbox.register("mutate", get_mutate_func(limits))
    
    best = toolbox.individual()
    update_fitnesses([best])
    
    time_on_plato = 0
    gen = 0
    while True:
        applicant = toolbox.clone(best)
        applicant = toolbox.mutate(applicant)
        del applicant.fitness.values
        
        update_fitnesses([applicant])

        new_best = toolbox.select([best, applicant])[0]
        
        if log:
            print("Gen {}: best: ({}) {}".format(gen, best.fitness.values[0], best))
            print("  applicant: ({}) {}".format(applicant.fitness.values[0], applicant))
        if callback is not None:
            if callback(best.fitness.values[0]) == True:
                break
        
        if new_best == best:
            time_on_plato += 1
            if time_on_plato == max_plato_time:
                break
        best = new_best
        gen += 1
        
    return best, best.fitness.values[0]

def one_plus_one(evaluate_func, limits, max_plato_time=100, log=False, callback=None):
    best, fit = None, None

    one_plus_one.stop = False
    def inner_callback(fitness):
        res = callback(fitness)
        one_plus_one.stop |= True if res == True else False
        return res

    while not one_plus_one.stop:
        new_best, new_fit = one_plus_one_onecall(evaluate_func, limits, max_plato_time, log, inner_callback)
        if fit == None or new_fit > fit:
            best, fit = new_best, new_fit

    return best, fit

def genetic(evaluate_func, limits, generations=300, popsize=20, elitepercent=.1, crossoveredpercent=.4, sbbx_eta=1.0, log=False, callback=None):
    assert elitepercent + crossoveredpercent <= 1.0
    elitesize = int(popsize * elitepercent)
    crossoveredsize = int(popsize * crossoveredpercent)
    crossoveredsize = crossoveredsize - 1 if crossoveredsize & 1 else crossoveredsize # make it even
    mutatedsize = popsize - elitesize - crossoveredsize
    assert mutatedsize >= 0

    toolbox.register("individual", tools.initIterate, creator.Individual, get_ind_function(limits))
    #toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    
    toolbox.register("evaluate", evaluate_func)
    toolbox.register("mutate", get_mutate_func(limits))
    toolbox.register("mate", lambda a, b: tools.cxSimulatedBinaryBounded(a, b, sbbx_eta, low=list(map(lambda t: t[0], limits)), up=list(map(lambda t: t[1], limits))))
    #toolbox.register("select", tools.selTournament, tournsize=3)
    #toolbox.register("select", tools.selBest)
    
    pop = [toolbox.individual() for _ in range(popsize)]
    
    #logbook = tools.Logbook()
    
    for gen in range(generations):
        update_fitnesses(pop)

        #record = stats.compile(pop)
        #logbook.record(gen=gen, **record)
        best=max(pop, key=lambda x: x.fitness.values[0])

        if log:
            print("Gen {}: best: {} !!! {}".format(gen, best, best.fitness.values[0]))#record['max']))
        if callback is not None:
            #if callback([(t, t.fitness.values[0]) for t in pop]) == True:
            if callback(best.fitness.values[0]) == True:
                break


        elite = tools.selBest(pop, k=elitesize)

        crossovered = list(map(toolbox.clone, tools.selTournament(pop, k=crossoveredsize, tournsize=2)))
        for child1, child2 in zip(crossovered[::2], crossovered[1::2]):
            toolbox.mate(child1, child2)
            del child1.fitness.values
            del child2.fitness.values
        
        mutated = list(map(toolbox.clone, random.sample(pop, mutatedsize)))
        for mutant in mutated:
            toolbox.mutate(mutant)
            del mutant.fitness.values
                
        pop[:] = elite + crossovered + mutated
        
    return best, best.fitness.values[0]#logbook