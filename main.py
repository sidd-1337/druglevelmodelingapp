import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

def calculate_drug_levels(initial_concentration, half_life_min, half_life_max, redose_schedule, duration):
    time_points = np.arange(0, duration + 1)
    concentrations_min = np.zeros_like(time_points, dtype=float)
    concentrations_max = np.zeros_like(time_points, dtype=float)
    concentrations_min[0] = initial_concentration
    concentrations_max[0] = initial_concentration

    redose_dict = {int(time): amount for time, amount in redose_schedule}

    for i in range(1, len(time_points)):
        decay_factor_min = 0.5 ** ((time_points[i] - time_points[i - 1]) / half_life_min)
        decay_factor_max = 0.5 ** ((time_points[i] - time_points[i - 1]) / half_life_max)

        concentrations_min[i] = concentrations_min[i - 1] * decay_factor_min
        concentrations_max[i] = concentrations_max[i - 1] * decay_factor_max

        if time_points[i] % 24 in redose_dict:
            concentrations_min[i] += redose_dict[time_points[i] % 24]
            concentrations_max[i] += redose_dict[time_points[i] % 24]

    return time_points, concentrations_min, concentrations_max

def find_local_extrema(concentrations):
    local_minima = argrelextrema(concentrations, np.less)[0]
    local_maxima = argrelextrema(concentrations, np.greater)[0]
    return local_minima, local_maxima

def main():
    st.title("Drug Level Modelling App")
    if 'redose_schedule' not in st.session_state:
        st.session_state.redose_schedule = []

    with st.sidebar:
        drug_name = st.text_input("Drug Name", value="Medication")
        concentration_unit = st.text_input("Concentration Unit (e.g., mg/L, mcg/mL)", value="mg/L")
        initial_concentration = st.number_input(f"Initial Concentration ({concentration_unit})", min_value=0.0,
                                                value=100.0)
        half_life_range = st.text_input("Half-Life Range (hours, e.g., 4.5-7)", value="6.0-12.0")

        try:
            half_life_min, half_life_max = map(float, half_life_range.split('-'))
            if half_life_min >= half_life_max:
                st.error("Minimum half-life must be less than maximum half-life.")
                return
        except ValueError:
            st.error("Invalid half-life range format. Please use the format 'min-max'.")
            return

        duration = st.number_input("Duration (hours)", min_value=1, value=72)

        st.header("Redosing Schedule")
        for i in range(len(st.session_state.redose_schedule)):
            cols = st.columns(2)
            with cols[0]:
                st.session_state.redose_schedule[i][0] = st.number_input(f"Hour {i + 1}",
                                                                         value=st.session_state.redose_schedule[i][0],
                                                                         key=f"hour_{i}")
            with cols[1]:
                st.session_state.redose_schedule[i][1] = st.number_input(f"Amount {i + 1} ({concentration_unit})",
                                                                         value=st.session_state.redose_schedule[i][1],
                                                                         key=f"amount_{i}")
            if st.button("Delete", key=f"del_{i}"):
                del st.session_state.redose_schedule[i]
                st.experimental_rerun()

        if st.button('Add Redose'):
            st.session_state.redose_schedule.append([0, 0])
            st.experimental_rerun()

    # Calculate drug levels whenever inputs change
    time_points, concentrations_min, concentrations_max = calculate_drug_levels(
        initial_concentration, half_life_min, half_life_max, st.session_state.redose_schedule, duration
    )

    # Find local minima and maxima
    local_minima_min, local_maxima_min = find_local_extrema(concentrations_min)
    local_minima_max, local_maxima_max = find_local_extrema(concentrations_max)

    fig, ax = plt.subplots()
    ax.plot(time_points, concentrations_min, label=f'{drug_name} Concentration (Min Half-Life)')
    ax.plot(time_points, concentrations_max, label=f'{drug_name} Concentration (Max Half-Life)')

    # Highlight local minima and maxima
    ax.scatter(time_points[local_minima_min], concentrations_min[local_minima_min], color='blue', label='Local Minima (Min Half-Life)')
    ax.scatter(time_points[local_maxima_min], concentrations_min[local_maxima_min], color='red', label='Local Maxima (Min Half-Life)')
    ax.scatter(time_points[local_minima_max], concentrations_max[local_minima_max], color='cyan', label='Local Minima (Max Half-Life)')
    ax.scatter(time_points[local_maxima_max], concentrations_max[local_maxima_max], color='magenta', label='Local Maxima (Max Half-Life)')

    ax.set_xlabel('Time (hours)')
    ax.set_ylabel(f'Concentration ({concentration_unit})')
    ax.set_title(f'{drug_name} Concentration Over Time')
    ax.legend()
    st.pyplot(fig)

if __name__ == "__main__":
    main()
