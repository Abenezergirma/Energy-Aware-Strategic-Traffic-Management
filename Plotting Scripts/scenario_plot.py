import folium
import json

# Central point of the original data (you can adjust this as needed)
original_center = [30, -90]

# New central point - DFW area
new_center = [32.8998, -97.0403]

# Calculate the difference
lat_diff = new_center[0] - original_center[0]
lon_diff = new_center[1] - original_center[1]

# Specify the name of the JSON file you want to read
filename = 'flight_plans.json'

# Read and load the JSON file

with open(filename, 'r') as file:
    flight_plans = json.load(file)

# Adjust the coordinates
for plan in flight_plans.values():
    for point in plan:
        point[0] += lat_diff
        point[1] += lon_diff

# Define the corner points and adjust them
corner_points = [
    [30.2, -90.1],  # Top left
    [29.8, -90.1],  # Bottom left
    [30, -89.6]     # Bottom right
]

adjusted_corners = [[point[0] + lat_diff, point[1] + lon_diff] for point in corner_points]

# Calculate the bounds of the rectangle
# Find the topmost, bottommost, leftmost, and rightmost points
topmost = max(adjusted_corners, key=lambda item: item[0])[0]
bottommost = min(adjusted_corners, key=lambda item: item[0])[0]
leftmost = min(adjusted_corners, key=lambda item: item[1])[1]
rightmost = max(adjusted_corners, key=lambda item: item[1])[1]

# Add a small offset to create space around the flight paths
offset_lat = 0.02  # Latitude offset
offset_lon = 0.03  # Longitude offset

# Define bounds using topmost, bottommost, leftmost, and rightmost with offsets
bounds = [
    [topmost + offset_lat, leftmost - offset_lon],  # Top left
    [bottommost - offset_lat, rightmost + offset_lon]  # Bottom right
]
        
# Create a map centered on the DFW area
# m = folium.Map(location=new_center, zoom_start=10)
m = folium.Map(location=new_center, zoom_start=10, tiles='CartoDB positron')
# m = folium.Map(location=new_center, zoom_start=10)

# Define new colors for different flight plans
colors = ['darkred', 'cadetblue', 'darkgreen']

# Add adjusted flight plans to the map
for idx, (key, coordinates) in enumerate(flight_plans.items()):
    # Adding departure and arrival icons
    # folium.Marker(coordinates[0], icon=folium.Icon(color=colors[idx], icon='plane', prefix='fa')).add_to(m)
    # folium.Marker(coordinates[-1], icon=folium.Icon(color=colors[idx], icon='flag', prefix='fa')).add_to(m)
    
    # Draw the flight path
    folium.PolyLine(
        locations=coordinates,
        color='blue',#cyan',#colors[idx],
        weight=2,
        opacity=1
    ).add_to(m)


# Add markers and circles for the intersection points
intersection_points = [[30, -90], [30, -89.8]]  # Two intersection points
for point in intersection_points:
    adjusted_point = [point[0] + lat_diff, point[1] + lon_diff]
    # folium.CircleMarker(adjusted_point, radius=4, color='black', fill=True, fill_color='black').add_to(m)
    # folium.Circle(
    #     adjusted_point,
    #     radius=600,  # radius in meters, adjust as needed
    #     color='cyan',
    #     fill=True,
    #     fill_opacity=0.0
    # ).add_to(m)

# Add an outer box (rectangle) with a contrasting color
# bounds = [[new_center[0]-0.2, new_center[1]-0.3], [new_center[0]+0.5, new_center[1]+0.3]]  # Adjust as needed
# folium.Rectangle(bounds, color='blue', fill=False, fill_opacity=0.3).add_to(m)

# Save to an HTML file
m.save('flight_plans_map.html')

