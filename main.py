import numpy as np
from tree import Tree
import random
import argparse
import copy
from struct import pack, unpack
from sys import getsizeof
from scipy.integrate import RK45, solve_ivp

data = Tree('data')
# data.show()  # отладка

def equetions(t, p):
    """
    t - моменты времени, но они не входят в уравнение
    p = P(L_i = l_i | T), (y в обозначениях из примеров)
    :return result: 1-D array (length = m * n)
    """
    result = []
    value_of_cur_fun = 0
    # по каждой популяции
    for pop in range(data.number_of_populations):
        # по всем образцам
        for i in range(data.samples_amount):
            # *** обращение к нужному p через индекс ( pop * data.samples_amount + i )
            ################################################
            # первое(1-я часть семмы) слагаемое
            for a in range(data.number_of_populations):
                value_of_cur_fun += data.migration_probability[a][pop] * p[a * data.samples_amount + i]
            ################################################
            # первое(2-я часть суммы) слагаемое
            for a in range(data.number_of_populations):
                value_of_cur_fun -= data.migration_probability[pop][a] * p[pop * data.samples_amount + i]
            ################################################
            # второе слагаемое
            first_multiplier = 0
            second_multiplier = 0
            sum_of_mult = 0
            for a in range(data.number_of_populations):
                first_multiplier = data.coalescence_probability[a] * p[a * data.samples_amount + i]
                for k in range(data.samples_amount):
                    if k != i:
                        second_multiplier += p[a * data.samples_amount + k]
                sum_of_mult += first_multiplier * second_multiplier
            value_of_cur_fun += p[pop * data.samples_amount + i] * sum_of_mult
            ################################################
            # третье слагаемое
            tmp_sum = 0
            for k in range(data.samples_amount):
                if k != i:
                    tmp_sum += p[pop * data.samples_amount + k]
            value_of_cur_fun -= p[pop * data.samples_amount + i] * data.coalescence_probability[pop] * tmp_sum
            ################################################
            result.append(value_of_cur_fun)
            value_of_cur_fun = 0

    return result


# по каждой популяции, и по всем образцам для каждой популяции
initial_states = data.get_initial_states()
# print(initial_states)  # отладка


# 0 - начальное время
# 15 - конечное время
for i in range(5):
    sol = RK45(equetions, t0=t, y0=initial_states, t_bound=15)
    initial_states = sol.f
    # print(sol.n)
    print(sol.t)
    print(sol.y)
    print(sol.f)
    # print(sol.direction)

"""
sol = solve_ivp(equetions, [0, 15], initial_states, t_eval=np.linspace(0, 15, 16))  # t_eval=np.linspace(0, 15, 16)
print(sol.t.reshape(-1,1))
print(sol.y[0])
print(sol.status)
"""
"""
sol = solve_ivp(lambda t, y: t-y, [0, 15], [2, 3],
                t_eval=np.linspace(0, 15, 16))
print(sol.t.reshape(-1,1))
print((sol.y))
"""
