# Dynamic-Movement-Primitives

Engineered a full assembly cycle for a SCARA robotic arm operating in a dynamic multi-conveyor environment. Computed
homogeneous transformations and forward kinematics to determine all intermediate end-effector poses. Planned
multi-segment trajectories in task space with continuous velocity profiles, accounting for scene geometry and
collision avoidance. Synchronized end-effector motion with two independently moving conveyors of different speeds
and phases. Extended the pipeline with Dynamic Movement Primitives (DMPs) for both
position and orientation (rotation about −z axis) to enable real-time adaptation to positional uncertainties (δx, δy)
and orientation perturbations (θδ) of target objects. Validated joint profiles, end-effector position/velocity
trajectories, including amplification near the folded posture. Kinematic simulation only, tested on three
perturbation cases using the MATLAB Robotics Toolbox.
