# Orbit Maneuver Calculator

This application calculates and visualizes different orbit maneuver types, including Hohmann transfers, Lambert transfers, and Bielliptic transfers. It also displays orbital parameters for initial, intermediate, and target orbits. The application is built using Streamlit and Plotly.
![alt text](https://github.com/JMMonte/Astrodynamics-visualizer-calculator/blob/main/Screenshot%202023-04-18%20at%2003.45.37.png?raw=true)
| | |
|:-------------------------:|:-------------------------:|
![alt text](https://github.com/JMMonte/Astrodynamics-visualizer-calculator/blob/main/Screenshot%202023-04-18%20at%2003.48.09.png?raw=true) | ![alt text](https://github.com/JMMonte/Astrodynamics-visualizer-calculator/blob/main/Screenshot%202023-04-18%20at%2003.47.47.png?raw=true) |


## Features

- Calculate and visualize Hohmann transfers, Lambert transfers, and Bielliptic transfers.
- Display orbital parameters for initial, intermediate, and target orbits.
- Visualize orbits in 3D using Plotly.
- Display ground track plots for Earth orbits.

## Requirements

- numpy
- plotly
- pandas
- streamlit
- astropy
- poliastro

## Installation

Clone the repository and install the required packages using pip:
Copy code

```bash
git clone https://github.com/your-username/orbit-maneuver-calculator.git cd orbit-maneuver-calculator pip install -r requirements.txt
```

To make sure you have all the dependencies installed, run:

```bash
pip install -r requirements.txt
```

## Usage

Run the Streamlit app:
Copy code

```bash
streamlit run app.py
```

Open your browser and go to <http://localhost:8501> to view the application.

## Application Structure

- app.py: Main file containing the Streamlit application.
- model.py: Contains functions and classes for orbit calculations and maneuver planning.

## License

This project is licensed under the MIT License.
