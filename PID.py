import math
class PID:
    def __init__(self, kp, ki, kd, dt):
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.dt = dt

        self.error = 0
        self.error_sum = 0
        self.error_diff = 0
        self.error_prev = 0

    def update(self, error):
        self.error = error
        self.error_sum += error * self.dt
        self.error_diff = (error - self.error_prev) / self.dt
        self.error_prev = error

        return self.kp * self.error + self.ki * self.error_sum + self.kd * self.error_diff
    
def find_cross_track_error(ref_path, state):
    # Find the nearest point on the path
    max_deviation = 3.5/2
    min_dist = 1000000
    min_index = 0
    for i in range(len(ref_path)):
        dist = (ref_path[i][0] - state[0])**2 + (ref_path[i][1] - state[1])**2
        if dist < min_dist:
            min_dist = dist
            min_index = i

    # Find the cross track error
    if min_index == 0:
        min_index += 1
    elif min_index == len(ref_path) - 1:
        min_index -= 1

    x1 = ref_path[min_index - 1][0]
    y1 = ref_path[min_index - 1][1]
    x2 = ref_path[min_index][0]
    y2 = ref_path[min_index][1]
    x3 = state[0]
    y3 = state[1]
    cross_track_error = (y2 - y1) * x3 - (x2 - x1) * y3 + x2 * y1 - y2 * x1
    cross_track_error /= math.sqrt((y2 - y1)**2 + (x2 - x1)**2)


    normalized_cross_track_error = cross_track_error / max_deviation

    return normalized_cross_track_error
