#%%

from JuLiFunc import geodistance
from SQL import get_clients_info, update_client_time, get_today_order_infos,get_sub_distance_mat
import pandas as pd
from TSP_alog import SiglemainDataFrame

def AssignLine(today_order_infos):
    """
    如果这个客户在数据库中，则将其分配到之前的线路，如果这个客户没在数据库中，则将其分配到里线路聚类中心最近的线路中去。
    :param credit_ids: 今日要配送的客户的信息
    :return: key=线路名字，value=线路所包含的客户，
    """
    # 从数据库中读取数据
    clients_data = get_clients_info()
    # key=线路名字，value=线路所包含的客户信息，客户信息的key=客户id,客户信息的value=客户经纬度
    LineInfos = {}
    for client_id in today_order_infos['order_client'].values:
        update_client_time(client_id)
        # 这个客户点的线路信息
        line_info = None
        # 如果数据库中有这个用户的信息，则将其分配到他之前所在的线路
        if client_id in clients_data['client_id'].values:
            line_info = clients_data[clients_data['client_id'] == client_id]['line_name'].values[0]
        else:
            # 取线路中心点
            center_credit = clients_data[clients_data['line_center'] == 1]
            min_distance = 1000000000
            # 当前客户点的经纬度
            cur_position = today_order_infos[today_order_infos['order_client'] == client_id]['order_position'].values[0]
            for center_id in center_credit['client_id'].values:
                # 里程距离
                cur_distance = geodistance(cur_position,
                                           clients_data[clients_data['client_id'] == center_id][
                                               'client_position'].values[
                                               0])
                if cur_distance < min_distance:
                    min_distance = cur_distance
                    line_info = clients_data[clients_data['client_id'] == center_id]['line_name'].values[0]
        # 获取该点的经纬度
        position = today_order_infos[today_order_infos['order_client'] == client_id]['order_position'].values[0]
        if line_info not in LineInfos:
            LineInfos[line_info] = {}
        # 将该点放入LineInfos中
        LineInfos[line_info][client_id] = position
    return LineInfos



def GetDistance(line_infos):
    for line_name,orders_info in line_infos.items():
        print(len(orders_info))
        sub_dis_mat = get_sub_distance_mat(orders_info,line_name)
        print(1)
        sub_dis_mat_index = pd.DataFrame(sub_dis_mat.index[1:])
        a,b,c=SiglemainDataFrame(sub_dis_mat,sub_dis_mat_index)
        print(c)
today_order_infos = get_today_order_infos()
LineInfos = AssignLine(today_order_infos)
GetDistance(LineInfos)
