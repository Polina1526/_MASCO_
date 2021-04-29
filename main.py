import functions as fun
import numpy as np
from tree import Tree
from scipy.integrate import solve_ivp


tree = Tree("coefficients", "migration_rates", "Newick_tree")
limits_list, lineage_list = fun.parse_tree_history(data=tree)

previous_states = None
sol_lines_list = []
sol_tree_list = []
for period in range(tree.samples_amount - 1):
    #########################################
    #                 ЛИНИИ                 #
    #########################################
    # тут находятся индексы образцов, принимающих участие в коалесценции в каждом периоде
    # при этом после кажой коалесценции инлексы после j-ого уменьшаются на один
    if period == 0:
        lineage = None
    else:
        lineage = lineage_list[period - 1]
    limits = limits_list[period]

    # последний столбец решений
    if period != 0:
        previous_states = sol_lines_list[-1].y[:, -1]
    initial_states = fun.create_initial(tree, period, previous_states=previous_states, lineage=lineage)  # ...

    # нужно ли находить значения функции только в точках с целым dt или в этом необходимости нет?
    t_span = np.linspace(limits[0], limits[1], 1001)  # limits[1] - limits[0] + 1)

    # решение диффуров для вероятностей линий на данном периоде
    sol_lines_list.append(solve_ivp(lambda t, y: fun.system_of_DE_for_lines(tree, y),
                                    t_span=limits,
                                    y0=initial_states,
                                    t_eval=t_span))

    # уменьшается cur_samples_amount для следующего периода после коалесцении
    tree.cur_samples_amount -= 1

print((fun.create_initial(tree, period=1, previous_states=sol_lines_list[-1].y[:, -1], lineage=lineage_list[-1]))[-1])
