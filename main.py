import numpy as np
from tree import Tree
import random
import argparse
import copy
from struct import pack, unpack
from sys import getsizeof
from scipy.integrate import solve_ivp

tree = Tree('data')
# отладка
# tree.show()


def equations(p: np.ndarray, data: Tree, period: int) -> np.ndarray:
    """
    :param p: p = P(L_i = l_i | T), (y в обозначениях из примеров)
    :param data:
    :return: result: 1-D array (length = m * n)
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
                second_multiplier = 0
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


"""
# по каждой популяции, и по всем образцам для каждой популяции
initial_states = tree.get_initial_states()
# отладка
# print(initial_states)
"""


def create_initial(data: Tree, period: int, previous_states: np.ndarray = None, lineage: np.ndarray = None) -> np.ndarray:
    """
    :param data:
    :param period:
    :param previous_states: probabilities before coalescence
    :param lineage: contains two samples that participated in coalescence
    :return: result: vector of initial states for current "period"
    """
    # *** обращение к нужному p через индекс ( pop * data.samples_amount + i )
    if period == 0:
        return data.get_initial_states()
    else:
        lineage = np.sort(lineage)
        result = []
        # cur_samples_sum = 0  ??нужно??

        # посчитаю все условные вероятности для lineage[0] и всех pop
        conditional_prob = []
        for pop in range(data.number_of_populations):
            conditional_prob.append(previous_states[pop * data.samples_amount + lineage[0]] *
                                    previous_states[pop * data.samples_amount + lineage[1]] *
                                    data.coalescence_probability[pop])
        sum_conditional_prob = np.sum(conditional_prob)  # ??? нужно ли переводить в массив ???
        # по каждой популяции
        for pop in range(data.number_of_populations):
            # по каждому образцу
            for i in range(data.samples_amount):
                # если образец участвовал в коалесценции и у него меньшее 'id', то есть он остался
                if i == lineage[0]:
                    # !!!!!!!!!!!!!! можно использовать уже посчитанные выше
                    result.append(previous_states[pop * data.samples_amount + i] *
                                  previous_states[pop * data.samples_amount + lineage[1]] *
                                  data.coalescence_probability[pop])
                # если образец участвовал в коалесценции и у него большее 'id', то есть он не остался
                elif i == lineage[1]:
                    continue
                # если образец не участвовал в коалесценции
                elif i != lineage[0] and i != lineage[1]:
                    result.append(previous_states[pop * data.samples_amount + i] * sum_conditional_prob)

            # cur_samples_sum += tree.number_of_samples[pop]  ??нужно??
        return result


def get_limits(period: int) -> np.ndarray:
    return


##########################################
# надо понять как это получать из дерева
limits_list = [np.array([0, 10]), np.array([10, 25])]
lineage_list = [np.ndarray([0, 1]), np.ndarray([0, 1])]
previous_states = None

sol_list = []
for period in range(tree.samples_amount - 1):

    lineage = lineage_list[period]  # ...
    initial_states = create_initial(tree, period, previous_states=previous_states, lineage=lineage)  # ...

    limits = limits_list[period]  # ...
    t_span = np.linspace(limits[0], limits[1], limits[1] - limits[0] + 1)
    sol_list.append(solve_ivp(lambda t, y: equations(y, tree, period), t_span=limits, y0=initial_states, t_eval=t_span))

    # уменьшается cur_samples_amount, так как следующий период после коалесцении
    tree.cur_samples_amount -= 1
    # последний столбец решений
    previous_states = sol_list[-1].y[:, -1]


for sol in sol_list:
    print(sol)
