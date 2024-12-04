# from pymodaq_femto.simulation import Simulator
#
# d = Simulator.params
# f = open("rst_params.rst", 'w')
# for param in d:
#     f.write('\n* **{}** '.format(param['title']))
#     if param['type'] == 'group':
#         for child in param['children']:
#             f.write('\n * {} '.format(child['title']))
#             f.write('``(type: {})`` '.format(child['type']))
#             f.write('Short description')
#     else:
#         f.write('(type: ``{}``) '.format(param['type']))
#         f.write('Short description\n')
#     # for k, v in value.items():
#     #     print(' * {}:'.format(k))
#     #     for i in v:
#     #         print('  * {}'.format(i))
#
# f.close()

from pymodaq_femto.retriever import Retriever

d = Retriever.params_in
f = open("ret_rst_params.rst", 'w')
for param in d:
    f.write('\n* **{}** '.format(param['title']))
    if param['type'] == 'group':
        for child in param['children']:
            f.write('\n * {} '.format(child['title']))
            f.write('``(type: {})`` '.format(child['type']))
            f.write('Short description')
            if child['type'] == 'group':
                for grandchild in child['children']:
                    f.write('\n  * {} '.format(grandchild['title']))
                    f.write('``(type: {})`` '.format(grandchild['type']))
                    f.write('Short description')
    else:
        f.write('(type: ``{}``) '.format(param['type']))
        f.write('Short description\n')
    # for k, v in value.items():
    #     print(' * {}:'.format(k))
    #     for i in v:
    #         print('  * {}'.format(i))

f.close()