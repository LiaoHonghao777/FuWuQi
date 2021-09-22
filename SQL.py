import datetime
import os

import pandas as pd
import pymysql
import time
from multiprocessing import Pool
from JuLiFunc import geodistance

db = pymysql.connect(host='10.16.127.139', port=3306, user='root', password='123456', db='wuliu')


def get_line_code():
    line_code = pd.read_sql("select * from line_code;", db)
    return line_code


def update_line_name():
    dirs = os.listdir("距离矩阵/")
    code = {}
    for dir in dirs:
        code1 = input("请输入英文代码{}.".format(dir))
        code[dir] = code1


def get_clients_info():
    clients_info = pd.read_sql("select * from client_infos;", db)
    return clients_info


def check_core_clients():
    """
    后台检查用户是否还是核心用户，如果当前时间离最后一次送货时间相隔超过30天，则取消他核心用户的标志
    """
    # 获取数据库中用户的id信息
    clients_ids = get_clients_info()['client_id']
    cursor = db.cursor()
    for client_id in clients_ids:
        # 获取当前的时间
        now_time_str = datetime.datetime.now().strftime("%Y-%m-%d")
        # 将当前时间转换为datetime
        now_time = datetime.datetime.strptime(now_time_str, "%Y-%m-%d")
        # 提取最后一次更新的时间
        pre_sql = "SELECT client_update_time from client_infos WHERE client_id='{}';".format(client_id)
        cursor.execute(pre_sql)
        pre_time = str(cursor.fetchone()[0])
        # 将最后一次更新时间转换为datetime
        pre_time = datetime.datetime.strptime(pre_time, "%Y-%m-%d")
        past_days = (now_time - pre_time).days
        if past_days > 30:
            sql = "UPDATE client_infos SET is_fre_client=0 WHERE client_id='{}';".format(now_time_str)
            cursor.execute(sql)
            db.commit()


def update_client_time(credit_id):
    """
    更新这位用户最后一次送货得时间。
    :param credit_ids:
    :return:
    """
    # 获取当前的时间
    cursor = db.cursor()
    now_time_str = datetime.datetime.now().strftime("%Y-%m-%d")
    sql = "UPDATE client_infos SET client_update_time='{}' WHERE client_id={};".format(now_time_str, credit_id)
    cursor.execute(sql)
    db.commit()


def get_distance_cangku():
    """
    获取客户id与仓库距离的数据表
    :return: cangku_distance---包含每个点和仓库的距离的表
    """
    cangku_distance = pd.read_sql("select * from cangku_distance;", db)
    return cangku_distance


def get_today_order_infos():
    """
    获取当日订单信息
    :return:
    """
    today_order_info = pd.read_sql("select * from today_order_info;", db)
    return today_order_info


def update_center_client(client_id):
    cursor = db.cursor()
    sql = "UPDATE client_infos SET line_center=1 WHERE client_id='{}';".format(client_id)
    cursor.execute(sql)
    db.commit()


def get_sub_distance_mat(line_credits_info, line_name):
    """
    依照当日客户信息,提取出当日所需要的距离矩阵
    :param line_credits_infos: 当日客户信息{客户id:客户经纬度}
    :param line_name: 线路名称
    :return: 返回当日所需要的距离矩阵
    """
    # 拼接数据库名字
    table_name = "{}_distance".format(line_name)
    distance_mat = pd.read_sql("select * from {};".format(table_name), db).values
    # 将cangku添加进客户信息中
    line_credits_info['cangku'] = '106.648824,29.468953'
    # 建立当日订单信息的距离矩阵
    sub_distance_mat = pd.DataFrame([])
    sub_distance_mat['cangku'] = -1
    sub_distance_mat.loc['cangku'] = -1
    for line_credit_id in line_credits_info.keys():
        if line_credit_id !='cangku':
            sub_distance_mat[str(line_credit_id)] = -1
            sub_distance_mat.loc[str(line_credit_id)] = -1
    for line in distance_mat:
        sub_distance_mat.loc[str(line[0]), str(line[1])] = line[2]
    for id1 in sub_distance_mat.columns:
        for id2 in sub_distance_mat.columns:
            if sub_distance_mat.loc[id1, id2] == -1:
                sub_distance_mat.loc[id1, id2] = geodistance(line_credits_info[id1], line_credits_info[id2])
    return sub_distance_mat


if __name__ == '__main__':
    pass
    # create_dis_table()
    # insert_dis_infos('距离矩阵/')
    # clients_info = get_clients_info()
    # print(clients_info.head())
    # check_core_clients()
    # create_today_order_info()
    # insert_today_order_info()
    # create_dis_table()
    # create_line_code()

    # distance = pd.read_sql("select * from {};".format("dj_distance"), db)
    # print(
    #     distance[(distance['start_position'] == 'cangku') & (distance['end_position'] == '8017561')]['distance'].values)
