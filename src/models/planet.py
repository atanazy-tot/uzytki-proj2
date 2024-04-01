class Planet:
    def __init__(self, name: str, mass: float, radius: float, rotation_speed: float):
        self.name = name
        self.mass = mass
        self.radius = radius
        self.rotation_speed = rotation_speed

    def __str__(self):
        return f"Planet {self.name}: Mass = {self.mass}, Radius = {self.radius}"