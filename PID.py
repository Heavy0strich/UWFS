
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