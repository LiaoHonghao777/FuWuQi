import pandas as pd
import random
from JuLiFunc import geodistance

def init_center(credit_ids: list, dictance_mat: pd.DataFrame, k: int):
    """
    Kmeans++的算法的初始点选择
    :param credit_ids:这条线路的客户id信息表
    :param dictance_mat: 这条线路的距离矩阵
    :param k: 聚类中心
    :return:初始聚类中心的id点
    """
    centers = []
    # 随机选择第一个点
    first_cneters = credit_ids[random.randint(0, len(credit_ids) - 1)]
    centers.append(first_cneters)
    while (len(centers) < k):
        # 每个点离已有聚类中心的距离
        dictance_x_to_center = {}
        for credit_id in credit_ids:
            if credit_id in centers:
                continue
            min_distance = 1000000
            # 寻找与每个中心点的最短距离
            for center in centers:
                distance = dictance_mat.loc[str(credit_id), str(center)]
                if distance < min_distance:
                    min_distance = distance
            dictance_x_to_center[credit_id] = min_distance
        # 计算每个点被选为下个中心点的概率
        sum_distance = 0
        for credit_id, distance in dictance_x_to_center.items():
            sum_distance += distance
        for credit_id in dictance_x_to_center.keys():
            dictance_x_to_center[credit_id] = dictance_x_to_center[credit_id] / sum_distance
        # 选取概率最高的作为初始中心点加入center
        max_p = -1
        max_p_credit_id = None
        for credit_id in dictance_x_to_center.keys():
            if dictance_x_to_center[credit_id] > max_p:
                max_p = dictance_x_to_center[credit_id]
                max_p_credit_id = credit_id
        centers.append(max_p_credit_id)
    return centers


def kmeanspp(credit_ids: list,id_position:dict, dictance_mat: pd.DataFrame,k: int):
    """
    通过kmeans++算法进行线路聚类中心的选择，其中迭代中心点修改为离质心最近的客户点
    :param credit_ids: 客户id信息表
    :param id_position: key为客户id,value为客户的经纬度
    :param dictance_mat: 距离矩阵
    :param k: 聚类数k
    :return: 返回代表该条线路的中心点
    """
    # 得到初始点
    centers = init_center(credit_ids, dictance_mat, k)
    # 迭代次数
    cnt = 0
    center_nodes = {}

    # 设置迭代次数为10次
    while cnt < 10:
        # 为每一个点划分所属类
        center_nodes = {}
        for credit_id in credit_ids:
            min_distance = 10000000000
            min_center = None
            for center in centers:
                distance = dictance_mat.loc[str(credit_id), str(center)]
                if distance < min_distance:
                    min_distance = distance
                    min_center = center
            if min_center not in center_nodes.keys():
                center_nodes[min_center] = [min_center]
            else:
                center_nodes[min_center].append(credit_id)
        # 进行迭代中心点的选择
        centers = []
        for center_id, center_sub_ids in center_nodes.items():
            center_id_x = 0
            center_id_y = 0
            for center_sub_id in center_sub_ids:
                center_sub_id_x = float(id_position[center_sub_id].split(',')[0])
                center_sub_id_y = float(id_position[center_sub_id].split(',')[1])
                center_id_x += center_sub_id_x
                center_id_y += center_sub_id_y
            center_id_x = center_id_x / len(center_sub_ids)
            center_id_y = center_id_y / len(center_sub_ids)
            cur_center_id = "{},{}".format(center_id_x,center_id_y)
            min_distance = 10000000000
            min_node = None
            for center_sub_id in center_sub_ids:
                distance = geodistance(cur_center_id,id_position[center_sub_id])
                if distance < min_distance:
                    min_distance = distance
                    min_node = center_sub_id
            centers.append(min_node)

        cnt = cnt + 1
    error = 0
    for center_id,center_ids in center_nodes.items():
        for client_id in center_ids:
            error+=geodistance(id_position[center],id_position[client_id])
    return centers,error

def cal_line_center(credit_ids: list,id_position:dict, distance_mat: pd.DataFrame):
    """
    计算每条线路的聚类中心,没有返回从数据库中读入
    :return:
    """
    min_error = 10000000000000
    min_center = None
    for k in range(3,10):
        center,error = kmeanspp(credit_ids,id_position,distance_mat,k)
        if error < min_error:
            min_error = error
            min_center = center
    return min_center
