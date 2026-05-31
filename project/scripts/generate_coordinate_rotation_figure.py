import numpy as np
import matplotlib.pyplot as plt


def get_rotation_matrix(alpha, beta, gamma):
    # Rotation angles in radians
    # R = Rz(alpha) * Ry(beta) * Rx(gamma)
    c1, s1 = np.cos(alpha), np.sin(alpha)
    c2, s2 = np.cos(beta), np.sin(beta)
    c3, s3 = np.cos(gamma), np.sin(gamma)

    Rz = np.array([[c1, -s1, 0], [s1, c1, 0], [0, 0, 1]])
    Ry = np.array([[c2, 0, s2], [0, 1, 0], [-s2, 0, c2]])
    Rx = np.array([[1, 0, 0], [0, c3, -s3], [0, s3, c3]])
    return Rz @ Ry @ Rx


def slerp(u, v, t, radius=0.6):
    u = u / np.linalg.norm(u)
    v = v / np.linalg.norm(v)
    dot = np.dot(u, v)
    theta = np.arccos(np.clip(dot, -1.0, 1.0))
    if np.abs(theta) < 1e-6:
        return radius * np.outer(1 - t, u)
    sin_theta = np.sin(theta)
    # Perform slerp for each t
    pts = []
    for ti in t:
        pt = (
            radius * (np.sin((1 - ti) * theta) * u + np.sin(ti * theta) * v) / sin_theta
        )
        pts.append(pt)
    return np.array(pts)


# Create figure with high DPI for premium publication quality
fig = plt.figure(figsize=(7, 7), dpi=300)
ax = fig.add_subplot(111, projection="3d")
ax.set_axis_off()

# Define orthonormal basis vectors for Catalog Frame F_C
origin = np.array([0, 0, 0])
x_C = np.array([1, 0, 0])
y_C = np.array([0, 1, 0])
z_C = np.array([0, 0, 1])

# Apply a mathematically exact SO(3) rotation matrix to guarantee orthonormal basis for F_S
# We choose distinct angles to make the 3D rotation visually clear and elegant
alpha = np.radians(35)  # 35 degrees around z
beta = np.radians(25)  # 25 degrees around y
gamma = np.radians(10)  # 10 degrees around x
R = get_rotation_matrix(alpha, beta, gamma)

x_S = R @ x_C
y_S = R @ y_C
z_S = R @ z_C

# Colors: sleek academic blue for Catalog, vibrant crimson for Surrogate
color_C = "#1f77b4"  # Catalog Frame (F_C)
color_S = "#d62728"  # Surrogate Frame (F_S)
color_arrow = "#ff7f0e"  # Rotation arrow/operator

# Draw Catalog Frame F_C axes (perfectly perpendicular basis)
ax.quiver(
    0,
    0,
    0,
    x_C[0],
    x_C[1],
    x_C[2],
    color=color_C,
    linewidth=2.5,
    arrow_length_ratio=0.08,
    pivot="tail",
)
ax.quiver(
    0,
    0,
    0,
    y_C[0],
    y_C[1],
    y_C[2],
    color=color_C,
    linewidth=2.5,
    arrow_length_ratio=0.08,
    pivot="tail",
)
ax.quiver(
    0,
    0,
    0,
    z_C[0],
    z_C[1],
    z_C[2],
    color=color_C,
    linewidth=2.5,
    arrow_length_ratio=0.08,
    pivot="tail",
)

# Draw Surrogate Frame F_S axes (guaranteed mathematically perpendicular because R is orthogonal)
ax.quiver(
    0,
    0,
    0,
    x_S[0],
    x_S[1],
    x_S[2],
    color=color_S,
    linewidth=2.5,
    arrow_length_ratio=0.08,
    pivot="tail",
)
ax.quiver(
    0,
    0,
    0,
    y_S[0],
    y_S[1],
    y_S[2],
    color=color_S,
    linewidth=2.5,
    arrow_length_ratio=0.08,
    pivot="tail",
)
ax.quiver(
    0,
    0,
    0,
    z_S[0],
    z_S[1],
    z_S[2],
    color=color_S,
    linewidth=2.5,
    arrow_length_ratio=0.08,
    pivot="tail",
)

# Offset for labels to avoid overlaps
offset = 0.08

# Labels for F_C axes
ax.text(
    *(x_C + offset * x_C),
    r"$\mathbf{x}_C$",
    color=color_C,
    fontsize=14,
    ha="center",
    va="center"
)
ax.text(
    *(y_C + offset * y_C),
    r"$\mathbf{y}_C$",
    color=color_C,
    fontsize=14,
    ha="center",
    va="center"
)
ax.text(
    *(z_C + offset * z_C),
    r"$\mathbf{z}_C$",
    color=color_C,
    fontsize=14,
    ha="center",
    va="center"
)

# Labels for F_S axes
ax.text(
    *(x_S + offset * x_S),
    r"$\mathbf{x}_S$",
    color=color_S,
    fontsize=14,
    ha="center",
    va="center"
)
ax.text(
    *(y_S + offset * y_S),
    r"$\mathbf{y}_S$",
    color=color_S,
    fontsize=14,
    ha="center",
    va="center"
)
ax.text(
    *(z_S + offset * z_S),
    r"$\mathbf{z}_S$",
    color=color_S,
    fontsize=14,
    ha="center",
    va="center"
)

# Label frame identifiers in the whitespace on the left half at the same height
ax.text(
    -0.3,
    0.15,
    0.9,
    r"Catalog Frame $\mathbf{F}_C$",
    color=color_C,
    fontsize=12,
    ha="center",
    va="center",
)
ax.text(
    -0.3,
    0.55,
    0.9,
    r"Surrogate Frame $\mathbf{F}_S$",
    color=color_S,
    fontsize=12,
    ha="center",
    va="center",
)

# Draw mathematically exact 3D rotation arc between z_C and z_S using SLERP
t_vals = np.linspace(0.0, 1.0, 100)
arc_pts = slerp(z_C, z_S, t_vals, radius=0.55)
ax.plot(
    arc_pts[:, 0],
    arc_pts[:, 1],
    arc_pts[:, 2],
    color=color_arrow,
    linestyle="-",
    linewidth=2,
)

# Add a clean arrowhead at the end of the arc
p_end = slerp(z_C, z_S, [0.90], radius=0.55)[0]
p_next = slerp(z_C, z_S, [1.0], radius=0.55)[0]
dp = p_next - p_end
dp = dp / np.linalg.norm(dp) * 0.05
ax.quiver(
    p_end[0],
    p_end[1],
    p_end[2],
    dp[0],
    dp[1],
    dp[2],
    color=color_arrow,
    arrow_length_ratio=0.4,
    pivot="tail",
    linewidth=2,
)

# Label the rotation operator R in the bottom-left whitespace
ax.text(
    -0.25,
    -0.25,
    0.0,
    r"$R \in \mathrm{SO}(3)$",
    color=color_arrow,
    fontsize=12,
    ha="center",
    va="center",
)

# Add small origin label
ax.text(-0.05, -0.05, -0.05, r"$O$", color="black", fontsize=12)

# Set equal aspect ratio to ensure no distortion (perpendicular axes remain visually perpendicular)
ax.set_xlim([-0.2, 1.1])
ax.set_ylim([-0.2, 1.1])
ax.set_zlim([-0.2, 1.1])
ax.set_box_aspect([1, 1, 1])

# Optimize 3D perspective angle
ax.view_init(elev=20, azim=30)

# Save to the project's figures folder with white background and tight margins
plt.savefig(
    "../figs/coordinate_frames_rotation.png",
    bbox_inches="tight",
    pad_inches=0.1,
    facecolor="white",
)
plt.close()
print("Figure generated successfully at project/figs/coordinate_frames_rotation.png")
