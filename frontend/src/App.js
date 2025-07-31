import React, { useState, useEffect } from 'react';
import { HashRouter as Router, Routes, Route } from 'react-router-dom';
import Dashboard from './pages/dashboard.jsx';
import CodeEditor from './pages/editor.jsx';
import Logs from './pages/logs.jsx';
import Sidebar from './components/Sidebar';
import './App.css';
import { ToastContainer } from 'react-toastify';
import 'react-toastify/dist/ReactToastify.css';

function App() {
  const [lightMode, setLightMode] = useState(true);
  const [matlabResults, setMatlabResults] = useState(null); // Results of MATLAB execution from CodeEditor
  const [message, setMessage] = useState(''); // Message to show status or error for CodeEditor

  // Toggle light mode on the body
  useEffect(() => {
    document.body.classList.toggle('light-mode', lightMode);
  }, [lightMode]);

  // handleRunMatlabCode is no longer needed here as a simulated function.
  // The editor.jsx will directly call setMatlabResults passed as a prop.

  return (
    <>
    <Router>
      <div className={`App dashboard-app ${lightMode ? 'light-mode' : ''}`}>
        <Sidebar lightMode={lightMode} setLightMode={setLightMode} />
        <div className="main">
          <Routes>
            <Route
              path="/"
              element={<Dashboard />}
            />
            <Route
              path="/dashboard"
              element={<Dashboard />}
            />
            <Route
              path="/editor"
              element={
                <CodeEditor
                  // No longer passing a simulated handleRunMatlabCode
                  setMatlabResults={setMatlabResults} // Still pass this to update App's state
                  setMessage={setMessage} // Pass this to allow editor to update the message
                />
              }
            />
            <Route
              path="/logs"
              element={<Logs />} />
          </Routes>

          {/* Display results from editor.jsx here */}
          <div style={{ padding: '20px', borderTop: '1px solid #ccc', marginTop: '20px' }}>
            <h3>Editor MATLAB Execution Output</h3>
            {message && <p>{message}</p>}
            {matlabResults ? (
              <pre style={{ backgroundColor: '#eee', padding: '10px', borderRadius: '5px', overflowX: 'auto' }}>
                {JSON.stringify(matlabResults, null, 2)}
              </pre>
            ) : (
              <p>No editor execution results yet.</p>
            )}
          </div>
        </div>
      </div>
    </Router>

    <ToastContainer position="top-right" autoClose={3000} />
    </>
  );
}

export default App;