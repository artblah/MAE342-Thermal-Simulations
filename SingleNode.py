import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Spacecraft thermal parameters
STEFAN_BOLTZMANN = 5.67e-8  # W/(m^2·K^4)
T_SPACE = 2.7  # Space temperature in Kelvin (cosmic background radiation)

class SurfaceComponent:
    """Represents a single surface component of the spacecraft"""
    def __init__(self, name, area, em, absBOL, absEOL, orientation):
        self.name = name 
        self.A = area #m^2
        self.epsilon = em
        self.alpha_BOL = absBOL
        self.alpha_EOL = absEOL
        self.orientation = orientation  # 'sun-facing', 'earth-facing', 'space-facing'
    
    def radiative_heat_loss(self, T):
        return self.epsilon * STEFAN_BOLTZMANN * self.A * (T**4 - T_SPACE**4)

class Louver:
    def __init__(self, area, emClosed, emOpen, Topen, Tclosed):
        self.A = area #m^2
        self.emClosed = emClosed
        self.emOpen = emOpen
        self.Topen = Topen
        self.Tclosed = Tclosed
        print(self.A, self.emClosed, self.emOpen, self.Topen, self.Tclosed)

    def effEmm(self, T):
        if T <= self.Tclosed:
            return self.emClosed
        elif T >= self.Topen:
            return self.emOpen
        return self.emOpen - (self.emOpen - self.emClosed) * (1 - T/self.Topen) / (1 - self.Tclosed/self.Topen)

    def radiative_heat_loss(self, T):
        return self.A * STEFAN_BOLTZMANN * (T**4 - T_SPACE**4) * self.effEmm(T)

class HeatInputs:
    def __init__(self, solar_flux, albedo_flux, earth_ir, internal_heat):
        self.solar_flux = solar_flux  # W/m^2
        self.albedo_flux = albedo_flux  # W/m^2
        self.earth_ir = earth_ir  # W/m^2
        self.internal_heat = internal_heat  # W

    def total_heat_input(self, components, BOL):
        total_heat = self.internal_heat
        for comp in components:
           if BOL:
                if comp.orientation == 'sun-facing':
                    total_heat += comp.alpha_BOL * self.solar_flux * comp.A
                elif comp.orientation == 'earth-facing':
                    total_heat += comp.alpha_BOL * self.albedo_flux * comp.A * 0.76
                    total_heat += comp.alpha_BOL * self.earth_ir * comp.A * 0.76
                elif comp.orientation == 'space-facing':
                    total_heat += comp.alpha_BOL * self.albedo_flux * comp.A * 0.2
                    total_heat += comp.alpha_BOL * self.earth_ir * comp.A * 0.2
           else:
                if comp.orientation == 'sun-facing':
                    total_heat += comp.alpha_EOL * self.solar_flux * comp.A
                elif comp.orientation == 'earth-facing':
                    total_heat += comp.alpha_EOL * self.albedo_flux * comp.A * 0.76
                    total_heat += comp.alpha_EOL * self.earth_ir * comp.A * 0.76
                elif comp.orientation == 'space-facing':
                    total_heat += comp.alpha_EOL * self.albedo_flux * comp.A * 0.2
                    total_heat += comp.alpha_EOL * self.earth_ir * comp.A * 0.2
        return total_heat

class SpacecraftNode:    
    def __init__(self, components, louver):
        self.components = components
        self.louver = louver
    
    def radiative_heat_loss(self, T):  #Kelvin
        total_loss = 0
        for component in self.components:
            total_loss += component.radiative_heat_loss(T)
        
        total_loss += self.louver.radiative_heat_loss(T)
        return total_loss
        
    
    def get_total_area(self):
        """Get total surface area of all components"""
        return sum(component.A for component in self.components)
    
    def print_summary(self):
        """Print summary of all components"""
        print(f"\n{'Component':<20} {'Area (m²)':<15} {'Emissivity':<15} {'Absorptivity BOL':<18} {'Absorptivity EOL':<18}")
        print("-" * 85)
        for comp in self.components:
            print(f"{comp.name:<20} {comp.A:<15.3f} {comp.epsilon:<15.3f} {comp.alpha_BOL:<18.3f} {comp.alpha_EOL:<18.3f}")
        
        # Print louver parameters
        print(f"{'Louver':<20} {self.louver.A:<15.3f} {self.louver.emClosed:<15.3f} {self.louver.emOpen:<18.3f}")
        print("-" * 85)
        print(f"{'TOTAL':<20} {self.get_total_area():<15.3f}")

    
    def equilibrium_temperature(self, heatInputs, BOL, Q):
        # Heat balance equation: Q_gen - Q_rad = 0
        heatInputs.internal_heat = Q  # Update internal heat generation for this calculation
        def heat_balance(T):
            return heatInputs.total_heat_input(self.components, BOL) - self.radiative_heat_loss(T)
        
        # Initial guess: room temperature
        T_initial_guess = 300
        T_eq = fsolve(heat_balance, T_initial_guess)[0]
        return T_eq


def main():
    # Define spacecraft surface components
    components = [
        SurfaceComponent("Solar Cells", area=0.258064, em=0.89, absBOL=0.112, absEOL = 0.322, orientation='sun-facing'),
        SurfaceComponent("Alum", area=0.129032, em=0.85, absBOL=0.9, absEOL = 0.945, orientation='space-facing'),
        SurfaceComponent("Antenna", area=0.00709676, em=0.88, absBOL=0.2, absEOL = 0.62, orientation='earth-facing'),
        SurfaceComponent("MLISun", area=0.629, em=0.01, absBOL=0.04, absEOL = 0.31, orientation='sun-facing'),
        SurfaceComponent("MLIEarth", area=0.856, em=0.01, absBOL=0.04, absEOL = 0.31, orientation='earth-facing'),
        SurfaceComponent("MLISpace", area=1.144, em=0.01, absBOL=0.04, absEOL = 0.31, orientation='space-facing')
    ]

    spacecraftLouver = Louver(area=0.462, emClosed=0.14, emOpen = 0.85, Topen = 40 + 273.15, Tclosed = -6 + 273.15)

    heatInputsPerigee = HeatInputs(solar_flux=1361, albedo_flux=63.5, earth_ir=245.5, internal_heat=0)
    heatInputsEclipse = HeatInputs(solar_flux=0, albedo_flux=0, earth_ir=0, internal_heat=0)  # Eclipse conditions

    # Create spacecraft node
    spacecraft = SpacecraftNode(components, spacecraftLouver)
  
    """
    # Print component summary
    print("Spacecraft Thermal Model - Component Summary:")
    spacecraft.print_summary()
    """
    

    # Range of heat generation values
    Q_gen_range = np.linspace(0, 200, 401)  # 0 to 200 W
    
    # Calculate equilibrium temperatures
    T_eq_range_BOL = np.zeros(Q_gen_range.size)
    for i, Q in enumerate(Q_gen_range):
        T_eq_range_BOL[i] = spacecraft.equilibrium_temperature(heatInputsPerigee, True, Q)

    T_eq_range_EOL = np.zeros(Q_gen_range.size)
    for i, Q in enumerate(Q_gen_range):
        T_eq_range_EOL[i] = spacecraft.equilibrium_temperature(heatInputsPerigee, False, Q)

    # Convert to Celsius for display
    T_eq_celsius_BOL = T_eq_range_BOL - 273.15
    T_eq_celsius_EOL = T_eq_range_EOL - 273.15

    # Create plot
    plt.style.use('classic')
    plt.figure(figsize=(10, 6))
    plt.plot(Q_gen_range, T_eq_celsius_BOL, 'b-', linewidth=2, label='BOL')
    plt.plot(Q_gen_range, T_eq_celsius_EOL, 'r-', linewidth=2, label='EOL')
    plt.axhline(y=-89, color='g', linestyle='--', linewidth=2)
    plt.axhline(y=89, color='g', linestyle='--', linewidth=2)
    plt.axhspan(-89, 89, alpha=0.2, color='green')
    plt.grid(True, alpha=0.3)
    plt.xlabel('Heat Generation (W)', fontsize=12)
    plt.ylabel('Equilibrium Temperature (°C)', fontsize=12)
    plt.title('Perigee Spacecraft Single Node Thermal Response', fontsize=12)
    plt.xlim(0, 200)
    plt.ylim(-100, 150)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.show()
    
    # Print some key values
    """print(f"\nSpacecraft Temperature at Different Heat Generations:")
    print(f"{'Heat Gen (W)':<15} {'T (K)':<15} {'T (°C)':<15}")
    print("-" * 45)
    for Q in [0, 20, 40, 60, 80,100, 120, 140]:
        T_k = spacecraft.equilibrium_temperature(heatInputsEclipse, BOL, Q)
        T_c = T_k - 273.15
        print(f"{Q:<15.0f} {T_k:<15.2f} {T_c:<15.2f}")
"""

if __name__ == "__main__":
    main()