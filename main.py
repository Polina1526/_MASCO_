import functions as fun
import numpy as np
from tree import Tree
import random
import argparse
import copy
from struct import pack, unpack
from sys import getsizeof
from scipy.integrate import solve_ivp

tree = Tree("coefficients", "migration_rates", "Newick_tree")
# отладка
# tree.show()

##############################################
# надо понять как это получать это из дерева #
##############################################
# Вариант 1, пример с рисунка из статьи
# limits_list = [np.array([0, 10]), np.array([10, 25])]
# lineage_list = [np.array([0, 1]), np.array([0, 1])]
# Вариант 2, пример из головы :)
limits_list = [np.array([0, 10]), np.array([10, 25]), np.array([25, 40])]
lineage_list = [np.array([1, 2]), np.array([0, 2]), np.array([0, 1])]


previous_states = None
sol_lines_list = []
sol_tree_list = []
for period in range(tree.samples_amount - 1):
    #########################################
    #                 ЛИНИИ                 #
    #########################################
    # тут находятся индексы образцов, принимающих участие в коалесценции в каждом периоде
    # при этом после кажой коалесценции инлексы после j-ого уменьшаются на один
    lineage = lineage_list[period]
    limits = limits_list[period]  # ...
    # последний столбец решений
    if period != 0:
        previous_states = sol_lines_list[-1].y[:, -1]
    initial_states = fun.create_initial(tree, period, previous_states=previous_states, lineage=lineage)  # ...

    # нужно ли находить значения функции только в точках с целым dt или в этом необходимости нет?
    t_span = np.linspace(limits[0], limits[1], limits[1] - limits[0] + 1)

    # решение диффуров для вероятностей линий на данном периоде
    sol_lines_list.append(solve_ivp(lambda t, y: fun.system_of_DE_for_lines(tree, y),
                                    t_span=limits,
                                    y0=initial_states,
                                    t_eval=t_span))

    ##########################################
    #                 ДЕРЕВО                 #
    ##########################################
    if period == 0:
        p_tree_before = None
    else:
        p_tree_before = sol_tree_list[-1].y[0][-1]
    tree_initial_state = fun.create_init_state_for_tree(tree, period,
                                                    p_tree_before=p_tree_before,
                                                    lines_prob=sol_lines_list[-1].y[:, -1],
                                                    lineage=lineage)

    # решение диффуров для вероятностей линий на данном периоде
    sol_tree_list.append(solve_ivp(lambda t, y: fun.DE_for_tree(tree, y, sol_lines_list[-1].y, t, limits[0]),
                                   t_span=limits,
                                   y0=tree_initial_state,
                                   t_eval=t_span))

    # уменьшается cur_samples_amount для следующего периода после коалесцении
    tree.cur_samples_amount -= 1


# for sol in sol_lines_list:
#     print(sol.y)

# for sol in sol_tree_list:
#     print(sol.y)

print(sol_tree_list[-1].y[0][-1])
