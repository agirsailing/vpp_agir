function [HS] = CalculateHydrostatics(ship,state)    
    [WS] = WetSections(ship, state);
    
    x = ship.x(1, :); %Vector transformation
   
    HS.V = trapz(x, (WS.Aws)); % Displacement volume
    HS.IWAx = trapz(x, (WS.bwl.^3) / 12); % Water plane area moment of inertia
    
    % Center of buoyancy coordinates
    HS.CoB(1) = trapz(x, x .* WS.Aws)/HS.V;
    HS.CoB(2) = trapz(x, WS.yB .* WS.Aws)/HS.V;
    HS.CoB(3) = trapz(x, WS.zB .* WS.Aws)/HS.V;  
    
    % Small heel angle stability assessment
    HS.BM0 = HS.IWAx / HS.V; 
    HS.KB0 = HS.CoB(3) + ship.KG; 
    %HS.KG0 = ship.CoG(3); 
    HS.KG0 = ship.KG;
    HS.GM0 = HS.KB0 + HS.BM0 - HS.KG0;
    
    % Hydrostatic lift force
    HS.FB0(1)=0; HS.FB0(2)=0;
    HS.FB0(3)=state.rau*state.g*HS.V; 

    % Transforming the coordinates of the CoB in the ship fixed coordinate system
    [HS.CoB0(1),HS.CoB0(2),HS.CoB0(3)] = TransShipfixedGlobal(HS.CoB(1),HS.CoB(2),HS.CoB(3),state.eta);

end