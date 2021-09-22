from math import radians, sin, cos, asin, sqrt


def geodistance(position1, position2):
    """
    计算两个经纬度直接的地图距离
    :param position1: 经纬度
    :param position2: 经纬度
    :return: 两个经纬度之间的地图距离
    """
    lng1, lat1 = position1.split(',')
    lng2, lat2 = position2.split(',')
    lng1, lat1, lng2, lat2 = map(radians, [float(lng1), float(lat1), float(lng2), float(lat2)])  # 经纬度转换成弧度
    dlon = lng2 - lng1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    distance = 2 * asin(sqrt(a)) * 6371 * 1000
    distance = round(distance, 3)
    return distance
