import React, { useState } from 'react';
import { Link, useLocation } from 'react-router-dom';

const Sidebar = ({ lightMode, setLightMode }) => {
  const location = useLocation();
  
  return (
    <div className="sidebar">
      <img src={`${process.env.PUBLIC_URL}/images/logo.png`} alt="Institute Logo" className="logo" />
      <nav>
        <ul>
          <li className={location.pathname === '/dashboard' ? 'active' : ''}>
            <Link to="/dashboard">
              <img src={`${process.env.PUBLIC_URL}/images/dashboard-light.png`} className="icon" alt="Dashboard Icon" /> 
              <span>Dashboard</span>
            </Link>
          </li>
          <li className={location.pathname === '/logs' ? 'active' : ''}>
            <Link to="/logs">
              <img src={`${process.env.PUBLIC_URL}/images/logs-light.png`} className="icon" alt="Logs Icon" /> 
              <span>Data & Logs</span>
            </Link>
          </li>
          <li className={location.pathname === '/editor' ? 'active' : ''}>
            <Link to="/editor">
              <img src={`${process.env.PUBLIC_URL}/images/code-light.png`} className="icon" alt="Code Editor Icon" /> 
              <span>Code Editor</span>
            </Link>
          </li>
        </ul>
      </nav>
      <div className="mode-toggle">
        <span>Light Mode</span>
        <label className="switch">
          <input 
            type="checkbox" 
            checked={lightMode}
            onChange={() => setLightMode(!lightMode)}
          />
          <span className="slider round"></span>
        </label>
      </div>
    </div>
  );
};

export default Sidebar; 
