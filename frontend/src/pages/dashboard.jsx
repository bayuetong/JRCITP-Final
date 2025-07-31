// src/pages/dashboard.jsx
import React, { useState, useEffect, useRef ,useCallback } from 'react';
import Plot from 'react-plotly.js'; // For rendering graphs
import Plotly from 'plotly.js-dist-min'; // For snapshot functionality
import { socket } from '../socket'; // Import the shared socket instance
import "../App.css"; // Your existing CSS
import { toast } from 'react-toastify';



const dropdownOptions = {
    'Radar Sensing': [
        { value: '', text: 'Select Graph' },
        { value: 'raw', text: 'Raw Breath' },
        { value: 'resp', text: 'Respiration' },
        { value: 'breathPSD', text: 'Breath PSD' },
        { value: 'heartbeatPSD', text: 'Heartbeat PSD' },
    ],
    'Vital Sign Detection': [
        { value: '', text: 'Select Graph' },
        { value: 'rangeDoppler', text: 'Range-Doppler' },
        { value: 'signalStrength', text: 'Signal Strength (N/A)' }, // Placeholder for now
        { value: 'clutterMap', text: 'Clutter Map (N/A)' },        // Placeholder for now
    ],
    'Human Vs Non-Human': [
        { value: '', text: 'Select Graph' },
        { value: 'humanVsNonHuman', text: 'Human vs Non-Human (N/A)' }, // Placeholder for now
    ],
};

const Dashboard = () => {
    // UI state
    const [selectedCategory, setSelectedCategory] = useState('Radar Sensing');
    const [graphOptions, setGraphOptions] = useState(dropdownOptions['Radar Sensing']);
    const [selectedGraph, setSelectedGraph] = useState('');
    const [scanStatus, setScanStatus] = useState('Inactive');
    const [uptime, setUptime] = useState('00:00:00');
    const [scanning, setScanning] = useState(false);
    const [paused, setPaused] = useState(false);

    // Mapping from dropdown 'value' to the 'name' property of the graph objects returned by MATLAB
const graphNameMapping = {
    'raw': 'JRC Raw Breath',
    'resp': 'JRC Filtered Breath',
    'breathPSD': 'JRC Breath PSD',
    'heartbeatPSD': 'JRC Heartbeat PSD',
    'rangeDoppler': 'JRC RDM',
    'sensorRaw': 'Sensor Raw Data',
    'sensorBreathPSD': 'Sensor Breath PSD',
    // 'signalStrength', 'clutterMap', 'humanVsNonHuman' are placeholders.
    // They won't display real data unless your MATLAB script generates them.
};

const createPlotlyTrace = useCallback((graphData, selectedGraphKey) => {
        console.log(`--- createPlotlyTrace Debug START ---`);
        console.log(`selectedGraphKey (from dropdown): ${selectedGraphKey}`);
        console.log(`graphData (full object from JSON):`, graphData);
        console.log(`graphData.name (from JSON): ${graphData?.name}`); // This is the display name

        if (!graphData || !graphData.name) {
            console.warn("createPlotlyTrace: Invalid graphData or missing 'name' property.");
            return null;
        }

        let trace = {};
        let layout = {
            title: { text: graphData.name || 'Graph Title' },
            autosize: true,
            margin: { t: 40, l: 60, r: 40, b: 60 },
            hovermode: 'closest',
            font: {
                size: 12
            }
        };

        switch (selectedGraphKey) {
            case 'rangeDoppler': // JRC RDM
                console.log('createPlotlyTrace: Detected Range-Doppler map.');
                if (graphData.rdm_data && graphData.rdm_data.z && graphData.rdm_data.x_axis && graphData.rdm_data.y_axis) {
                    const z = graphData.rdm_data.z;
                    const x = graphData.rdm_data.x_axis;
                    const y = graphData.rdm_data.y_axis;

                    if (!z || z.length === 0 || !z[0] || z[0].length === 0) {
                        console.error('createPlotlyTrace: RDM z data is empty or invalid.');
                        return null;
                    }

                    trace = {
                        z: z,
                        x: x,
                        y: y,
                        type: 'heatmap',
                        colorscale: 'Jet',
                        colorbar: { title: 'Intensity' },
                        name: graphData.name,
                    };
                    layout.xaxis = { title: 'Velocity (m/s)', autorange: true };
                    layout.yaxis = { title: 'Range (m)', autorange: true };
                    console.log('createPlotlyTrace: Successfully created RDM trace and layout.');
                } else {
                    console.warn(`createPlotlyTrace: Missing RDM data fields (z, x_axis, y_axis) for ${graphData.name}.`, graphData);
                    return null;
                }
                break;

            case 'raw': // JRC Raw Breath (time-series)
            case 'resp': // JRC Filtered Breath (time-series)
            case 'sensorRaw': // Sensor Raw Data (time-series)
                console.log(`createPlotlyTrace: Detected Time-series graph: ${graphData.name}`);
                if (graphData.time && graphData.signal) {
                    trace = {
                        x: graphData.time,
                        y: graphData.signal,
                        type: 'scatter',
                        mode: 'lines',
                        name: graphData.name,
                        line: { color: '#0074D9' },
                    };
                    layout.xaxis = { title: 'Time (s)', autorange: true };
                    layout.yaxis = { title: 'Amplitude', autorange: true };
                    console.log('createPlotlyTrace: Successfully created Time-series trace and layout.');
                } else {
                    console.warn(`createPlotlyTrace: Missing 'time' or 'signal' for ${graphData.name}.`, graphData);
                    return null;
                }
                break;

            case 'breathPSD': // JRC Breath PSD (frequency/PSD)
            case 'heartbeatPSD': // JRC Heartbeat PSD (frequency/PSD)
            case 'sensorBreathPSD': // Sensor Breath PSD (frequency/PSD)
                console.log(`createPlotlyTrace: Detected PSD graph: ${graphData.name}`);
                if (graphData.frequency && graphData.psd) {
                    trace = {
                        x: graphData.frequency,
                        y: graphData.psd,
                        type: 'scatter',
                        mode: 'lines',
                        name: graphData.name,
                        line: { color: '#FF4136' },
                    };
                    layout.xaxis = { title: 'Frequency (Hz)', autorange: true };
                    layout.yaxis = { title: 'Power/Magnitude', type: 'linear', autorange: true };
                    console.log('createPlotlyTrace: Successfully created PSD trace and layout.');
                } else {
                    console.warn(`createPlotlyTrace: Missing 'frequency' or 'psd' for ${graphData.name}.`, graphData);
                    return null;
                }
                break;

            default:
                console.warn(`createPlotlyTrace: Unsupported graph type or missing expected data for '${graphData.name}' (selected key: ${selectedGraphKey}).`, graphData);
                return null;
        }

        console.log(`--- createPlotlyTrace Debug END ---`);
        return { data: [trace], layout: layout };
    }, []);

    // Timer refs
    const startTimeRef = useRef(null);
    const pauseStartRef = useRef(null);
    const pausedDurationRef = useRef(0);
    const intervalRef = useRef(null);

    // State for MATLAB Analysis Results
    const [analysisGraphs, setAnalysisGraphs] = useState([]); // Array of Plotly-ready graph objects
    const [analysisResults, setAnalysisResults] = useState(null); // Scalar results (e.g., heart rate, breath rate)
    const [analysisError, setAnalysisError] = useState(null); // Error message from MATLAB
    const [isLoading, setIsLoading] = useState(false); // Loading indicator for analysis
    const [mainPlotlyData, setMainPlotlyData] = useState(null);

    // --- Graph Option/Category Management ---
    // Updates graph options when the selected category changes.
    useEffect(() => {
        setGraphOptions(dropdownOptions[selectedCategory]);
        setSelectedGraph(''); // Reset selected graph when category changes
    }, [selectedCategory]);

    // --- Uptime Timer Management ---
    // Manages the real-time uptime display.
    useEffect(() => {
        if (scanning && !paused) {
            intervalRef.current = setInterval(updateUptime, 1000);
        } else {
            clearInterval(intervalRef.current);
        }
        return () => clearInterval(intervalRef.current); // Cleanup on unmount or dependency change
    }, [scanning, paused]);

    // Calculates and updates the uptime string.
    const updateUptime = () => {
        if (!startTimeRef.current) return;
        const now = new Date();
        const elapsed = now - startTimeRef.current - pausedDurationRef.current;
        const hours = String(Math.floor(elapsed / 3600000)).padStart(2, '0');
        const minutes = String(Math.floor((elapsed % 3600000) / 60000)).padStart(2, '0');
        const seconds = String(Math.floor((elapsed % 60000) / 1000)).padStart(2, '0');
        setUptime(`${hours}:${minutes}:${seconds}`);
    };

    // --- Socket.IO Listener for Radar Analysis Results ---
    // Listens for 'radar_analysis_results' events from the backend.
    useEffect(() => {
        const handleRadarAnalysisResults = (data) => {
            setIsLoading(false); // Analysis has finished
            console.log('--- RECEIVED RADAR ANALYSIS RESULTS ---');
            console.log('Raw data received:', data);
            console.log('Type of data.graphs:', typeof data.graphs, data.graphs);

            if (data.error) {
                console.error('MATLAB Radar Analysis Error:', data.error);
                setAnalysisError(data.error);
                setAnalysisGraphs([]); // Clear graphs on error
                setAnalysisResults(null); // Clear results on error
                setMainPlotlyData(null); // Clear main plot on error
            } 
            
            else 
            
            {
                
                setAnalysisGraphs(data.graphs);
                console.log('analysisGraphs state set to:', data.graphs);

                const expectedGraphName = graphNameMapping[selectedGraph];
                console.log('Selected Graph from Dropdown:', selectedGraph);
                console.log('Expected Graph Name (from mapping):', expectedGraphName);
                
                const foundGraphData = data.graphs.find(g => g.name === expectedGraphName);
                console.log('Found graph data from results:', foundGraphData);

                if (foundGraphData) 
                {
                    const renderData = createPlotlyTrace(foundGraphData, selectedGraph);
                    console.log('Render data generated by createPlotlyTrace:', renderData);
                    setMainPlotlyData(renderData);
                    console.log('Main plot data state updated to:', renderData);
                }
                
                else 
                {
                    setMainPlotlyData(null);
                    console.warn(`No graph data found for selected graph: ${expectedGraphName}`);
                }
                setAnalysisResults(data.results);
                setAnalysisError(null);
            }
        }

        socket.on('radar_analysis_results', handleRadarAnalysisResults);

        return () => {
            socket.off('radar_analysis_results', handleRadarAnalysisResults);
        };
    }, [selectedGraph, createPlotlyTrace]);

	const loadRadarDataFromFile = useCallback(async () => {
    setIsLoading(true);
    setAnalysisError(null);
    setMainPlotlyData(null);

    try {
        const response = await fetch('http://localhost:5000/api/get_latest_radar_data');
        if (!response.ok) {
            const errorText = await response.text();
            throw new Error(`HTTP error! status: ${response.status}, message: ${errorText}`);
        }
        const data = await response.json();
        console.log('Data loaded from JSON file:', data);

        if (Array.isArray(data) && data.length > 0) {
            const allGraphData = data.map(entry => entry.graph_data).filter(Boolean);
            const firstEntry = data[0];
            const metadata = firstEntry.metadata || {};
            const results = firstEntry.results || {};

            setAnalysisGraphs(allGraphData);
            setAnalysisResults(results);

            const expectedGraphName = graphNameMapping[selectedGraph];
            const foundGraphData = allGraphData.find(g => g.name === expectedGraphName);

            if (foundGraphData) {
                const renderData = createPlotlyTrace(foundGraphData, selectedGraph);
                setMainPlotlyData(renderData);
                console.log(`Initial plot for ${selectedGraph} generated and set after file load.`);
            } else {
                setMainPlotlyData(null);
                console.warn(`No graph data found for selected graph: ${expectedGraphName} in the loaded file.`);
            }
        } else {
            
            setAnalysisError("JSON file does not contain an expected array of data objects, or is empty.");
        }

    } catch (error) {
        console.error('Error loading radar data from file:', error);
        setAnalysisError(`Failed to load data: ${error.message}`);
        toast.info('Error loading data!');
    } finally {
        setIsLoading(false);
    }
}, [selectedGraph, createPlotlyTrace]);


// Re-render the graph whenever the user selects a different graph from the dropdown
// or when new analysis data is received from the backend.
useEffect(() => {
    // Exit early if no graph is selected or if analysis data hasn't been received yet
    if (!selectedGraph || analysisGraphs.length === 0) return;

    // Get the actual graph name that matches the selected dropdown label
    const expectedGraphName = graphNameMapping[selectedGraph];

    // Search for the corresponding graph data in the list of received analysis graphs
    const foundGraphData = analysisGraphs.find(g => g.name === expectedGraphName);

    if (foundGraphData) {
        // Convert the found graph data into Plotly-compatible trace(s) and update the state
        const renderData = createPlotlyTrace(foundGraphData, selectedGraph);
        setMainPlotlyData(renderData);
    } else {
        // No matching graph found — clear the graph area and log a warning
        setMainPlotlyData(null);
        console.warn(`No graph found in analysisGraphs for selected graph: ${expectedGraphName}`);
    }
}, [selectedGraph, analysisGraphs, createPlotlyTrace]);


    // --- Function to initiate Radar Analysis ---
    // Emits a 'run_radar_analysis' event to the backend with configuration data.
    const runRadarAnalysis = async () => {
        setIsLoading(true);         // Indicate loading
        setAnalysisError(null);     // Clear previous errors
        setAnalysisGraphs([]);      // Clear previous graph data
        setAnalysisResults(null);   // Clear previous scalar results

        // Configuration data sent to MATLAB (these keys MUST match expected fields in runAnalysisWrapper.m)
        const configData = {
            runTests: true,
            sensor_test: true,
            hearbeat_test: true, // This appears to be a typo from the original, assuming it means heartbeat_test
            speed_test: false,
            real_time: true,
            ts: 0.06,
            // IMPORTANT: These paths must be ABSOLUTE paths on the SERVER (Flask) machine
            img_filename: '/home/xzhang/Desktop/matlab/images/jarcat.jpg',
            strCOM: '/dev/ttyACM1',
        };
        console.log('Emitting run_radar_analysis with config:', configData);
        socket.emit('run_radar_analysis', configData);
        // The results will be received by the 'radar_analysis_results' useEffect listener
    };

    // --- Control Button Handlers ---
    const handleStart = () => {
        if (!scanning) {
            startTimeRef.current = new Date();
            pausedDurationRef.current = 0;
            //runRadarAnalysis(); // Trigger the MATLAB analysis
            loadRadarDataFromFile();
            toast.info('Graph loading started!');
        } else if (paused) {
            pausedDurationRef.current += new Date() - pauseStartRef.current;
        }
        setScanning(true);
        setPaused(false);
        setScanStatus('Active');
    };

    const handlePause = () => {
        if (scanning && !paused) {
            pauseStartRef.current = new Date();
            setPaused(true);
            setScanStatus('Paused');
            toast.info('Graph loading paused!')
        }
    };

    const handleStop = () => {
        setScanning(false);
        setPaused(false);
        startTimeRef.current = null;
        pauseStartRef.current = null;
        pausedDurationRef.current = 0;
        setScanStatus('Stopped');
        setUptime('00:00:00');
        // Clear displayed analysis results when stopped
        setAnalysisGraphs([]);
        setAnalysisResults(null);
        setAnalysisError(null);
        setIsLoading(false);
        toast.info('Graph loading stopped!')
    };

    // --- Render Function ---
    return (
        <div className="content">
            <div className="left-main">
                <div className="config-panel">
                    <h3>Configuration</h3>
                    {Object.keys(dropdownOptions).map(cat => (
                        <div className="config-option" key={cat}>
                            <button
                                className={selectedCategory === cat ? 'active' : ''}
                                onClick={() => setSelectedCategory(cat)}>{cat}</button>
                        </div>
                    ))}
                    <div className="config-divider"></div>
                    <h3>System Status</h3>
                    <div className="config-status">
                        <div className="status-item"><span className="label">Scanning:</span> <span>{scanStatus}</span></div>
                        {/* Signal Strength will need to be populated from analysisResults if MATLAB provides it */}
                        <div className="status-item"><span className="label">Signal Strength:</span> <span>-</span></div>
                        <div className="status-item"><span className="label">Uptime:</span> <span>{uptime}</span></div>
                    </div>
                    {/* Display Scalar Analysis Results */}
                    <div className="config-divider"></div>
                    <h3>Analysis Results</h3>
                    {isLoading && <p>Running analysis...</p>}
                    {analysisError && <p style={{ color: 'red' }}>Error: {analysisError}</p>}
                    
                    {analysisResults ? (
                        <div className="analysis-results">
                            <p><strong>JRC Breath Rate:</strong>{' '}{analysisResults.jrcBreathRate != null ? analysisResults.jrcBreathRate.toFixed(2): 'N/A'}{' '}times/min</p>
                            {/* Check if heart rate exists before displaying */}
                            {analysisResults.jrcHeartRate != null && ( <p><strong>JRC Heart Rate:</strong> {analysisResults.jrcHeartRate.toFixed(2)}{' '}times/min</p>
                            )}
                            {analysisResults.sensorBreathRate !== undefined && (
                                <p><strong>Sensor Breath Rate:</strong>{' '}{analysisResults.sensorBreathRate != null ? analysisResults.sensorBreathRate.toFixed(2): 'N/A'}{' '}times/min</p>
                            )}
                            <p><strong>Range:</strong> {' '}{analysisResults.range != null ? analysisResults.range.toFixed(2) : 'N/A'}{' '} m</p>
                            <p><strong>Speed:</strong> {' '}{analysisResults.speed != null ? analysisResults.speed.toFixed(2): 'N/A'}{' '}m/s</p>
                        </div>
                    ) : (
                        // Show message only if not loading
                        !isLoading && <p>No analysis results yet.</p>
                    )}
                </div>
            </div>
            <div className="right-main">
                <div className="dashboard-row">
                    <div className="status-panel">
                        <div className="button-group">
                            <button className="control-button" onClick={handleStart} disabled={isLoading}>
                                <img src={`${process.env.PUBLIC_URL}/images/start.png`} className="icon" alt="Start" />
                                <span>{isLoading ? 'Starting...' : 'Start'}</span>
                            </button>
                            <button className="control-button" onClick={handlePause} disabled={!scanning || paused}>
                                <img src={`${process.env.PUBLIC_URL}/images/pause.png`} className="icon" alt="Pause" />
                                <span>Pause</span>
                            </button>
                            <button className="control-button" onClick={handleStop} disabled={!scanning && !isLoading}>
                                <img src={`${process.env.PUBLIC_URL}/images/stop.png`} className="icon" alt="Stop" />
                                <span>Stop</span>
                            </button>
                        </div>
                        <div className="status-label">Program Status</div>
                    </div>
                    <section className="graph-panel">
                        <h2>Graph Selector</h2>
                        <select value={selectedGraph} onChange={e => setSelectedGraph(e.target.value)} disabled={isLoading}>
                            {graphOptions.map(opt => (
                                <option key={opt.value} value={opt.value}>{opt.text}</option>
                            ))}
                        </select>
                    </section>
                </div>
                <section className="display-panel">
                    <div className="top-row">
                        <h2 id="graphTitle">
                            {/* Display the title from the generated Plotly layout */}
                            {mainPlotlyData?.layout?.title?.text || 'Graph'}
                        </h2>
                        <div className="action-buttons">
                            <button
                                className="snapshot"
                                onClick={async () => {
                                    const graphArea = document.getElementById("graphArea");
                                    // Ensure a plot element and data exist to snapshot
                                    if (!graphArea || !mainPlotlyData || !mainPlotlyData.trace) {
                                        alert("No graph to snapshot. Please select a graph and run analysis first.");
                                        return;
                                    }
                                    try {
                                        // Use Plotly.toImage to capture the graph as a PNG image
                                        const imageData = await Plotly.toImage(graphArea, { format: "png", scale: 2 }); // Scale for higher resolution
                                        const logs = JSON.parse(localStorage.getItem("logs") || "[]");
                                        logs.push({
                                            timestamp: new Date().toISOString().slice(0, 19).replace("T", " "),
                                            type: "SNAPSHOT",
                                            message: `Snapshot captured for ${mainPlotlyData.layout.title.text}`,
                                            data: {
                                                graphType: mainPlotlyData.layout.title.text,
                                                uptime,
                                                scanStatus,
                                                imageData // Base64 encoded image data of the plot
                                            }
                                        });
                                        localStorage.setItem("logs", JSON.stringify(logs));
                                        alert("Snapshot saved!");
                                    } catch (e) {
                                        console.error("Error capturing snapshot:", e);
                                        alert("Failed to capture snapshot. Ensure a graph is displayed and loaded.");
                                    }
                                }}
                                disabled={!mainPlotlyData} // Disable button if no graph data is loaded
                            >
                                <img src={`${process.env.PUBLIC_URL}/images/Camera-light.png`} className="icon" alt="Camera" />
                                Snapshot
                            </button>
                            <button
                                className="snapshot"
                                onClick={() => {
                                    // Task 2: "All Graphs" Button Functionality
                                    // 1. Store analysisGraphs in localStorage for access by the new window
                                    if (analysisGraphs && analysisGraphs.length > 0) {
                                        try {
                                            localStorage.setItem('allAnalysisGraphs', JSON.stringify(analysisGraphs));
                                            console.log('Analysis graphs stored in localStorage for new window.');
                                        } catch (e) {
                                            console.error('Failed to store analysisGraphs in localStorage:', e);
                                            alert('Failed to prepare data for all graphs. Please try again.');
                                            return;
                                        }
                                    } else {
                                        alert('No analysis results available to display all graphs. Please run an analysis first.');
                                        return;
                                    }

                                    // Log the action to localStorage
                                    const logs = JSON.parse(localStorage.getItem("logs") || "[]");
                                    logs.push({
                                        timestamp: new Date().toISOString().slice(0, 19).replace("T", " "),
                                        type: "LOG",
                                        message: `Opened "All Graphs" window`,
                                        data: {
                                            uptime,
                                            scanStatus,
                                        }
                                    });
                                    localStorage.setItem("logs", JSON.stringify(logs));

                                    // Open a new, blank window
                                    const graphWindow = window.open('', '_blank', 'width=1200,height=800'); // Increased size for better display
                                    if (!graphWindow) {
                                        alert("Please allow pop-ups for this site to view all graphs.");
                                        return;
                                    }
                                    // Write the HTML content directly into the new window
                                    graphWindow.document.write(`
                                        <html>
                                            <head>
                                                <title>All Radar Analysis Graphs</title>
                                                <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
                                                <style>
                                                    body { font-family: Arial, sans-serif; margin: 0; padding: 20px; background-color: #f0f2f5; }
                                                    h1 { text-align: center; margin-bottom: 20px; color: #333; }
                                                    #plots {
                                                        display: grid;
                                                        grid-template-columns: repeat(auto-fit, minmax(450px, 1fr)); /* Responsive grid */
                                                        gap: 30px; /* Space between plots */
                                                        justify-content: center;
                                                        padding: 0 10px; /* Some padding on the sides */
                                                    }
                                                    .plot-container {
                                                        background-color: #fff;
                                                        border: 1px solid #ddd;
                                                        padding: 15px;
                                                        box-shadow: 0 4px 12px rgba(0,0,0,0.1);
                                                        border-radius: 10px;
                                                        display: flex;
                                                        flex-direction: column;
                                                        align-items: center;
                                                        justify-content: center;
                                                        min-height: 350px; /* Minimum height for each plot */
                                                        width: 100%; /* Ensures it fits within grid column */
                                                        box-sizing: border-box; /* Include padding and border in the element's total width and height */
                                                    }
                                                    /* Plotly's internal div. This is crucial for responsive plots within containers */
                                                    .plot-container .js-plotly-plot { width: 100% !important; height: 100% !important; }
                                                    .no-data-message {
                                                        text-align: center;
                                                        color: #666;
                                                        font-style: italic;
                                                        margin-top: 50px;
                                                    }
                                                </style>
                                            </head>
                                            <body>
                                                <h1>All Radar Analysis Graphs</h1>
                                                <div id="plots"></div>
                                                <script>
                                                    // Mapping from MATLAB's graph name to the type string expected by createPlotlyTraceForNewWindow
                                                    const graphNameMapping = {
                                                        'raw': 'JRC Raw Breath',
                    					                'resp': 'JRC Filtered Breath',
                    					                'breathPSD': 'JRC Breath PSD',
                    					                'heartbeatPSD': 'JRC Heartbeat PSD',
                    					                'rangeDoppler': 'JRC RDM',
                    					                'sensorRaw': 'Sensor Raw Data',
                    					                'sensorBreathPSD': 'Sensor Breath PSD',
                                                    };
							 // Helper to flatten nested arrays (if MATLAB sends [value] instead of value)
                					const flattenArray = (arr) => Array.isArray(arr) ? arr.map(val => (Array.isArray(val) && val.length > 0) ? val[0] : val) : [];

                                                    // This function is a direct copy of createPlotlyTrace from Dashboard.jsx
                                                    // It needs to be self-contained in the new window.
                                                    function createPlotlyTraceForNewWindow(graphData) {
                                                        if (!graphData || !graphData.name) {
                                                            console.warn("Invalid graphData provided to createPlotlyTraceForNewWindow.");
                                                            return null;
                                                        }

							let trace = { name: graphData.name };
                    					let layout = {
                        				title: { text: graphData.name || 'Graph Title', font: { size: 16 } },
                        				autosize: true,
                        				margin: { t: 50, l: 70, r: 40, b: 70 }, // Adjusted margins for better fit
                        				hovermode: 'closest',
                        				paper_bgcolor: '#f8f9fa',
                        				plot_bgcolor: '#ffffff',
                        				font: { color: '#333' }
                    					};
                                                        // Determine graphType based on the graphData.name provided by MATLAB
                                                        // Try direct lookup, or try removing "JRC " prefix
                                                        const internalGraphTypeKey = Object.keys(graphNameMapping).find(key => graphNameMapping[key] === graphData.name);

if (internalGraphTypeKey === 'rangeDoppler' && graphData.rdm_data?.z) {
                        const z_data = graphData.rdm_data.z;
                        const x_data = flattenArray(graphData.rdm_data.x_axis);
                        const y_data = flattenArray(graphData.rdm_data.y_axis);

                        if (!z_data || z_data.length === 0 || !z_data[0] || z_data[0].length === 0) {
                            console.error('createPlotlyTraceForNewWindow: RDM z data is empty or invalid.');
                            return null;
                        }

                        trace = {
                            ...trace,
                            z: z_data,
                            x: x_data,
                            y: y_data,
                            type: 'heatmap',
                            colorscale: 'Viridis',
                            colorbar: { title: 'Intensity (dB)', titleside: 'right' },
                            hovertemplate: 'Range: %{y:.2f} m<br>Doppler: %{x:.2f} m/s<br>Intensity: %{z:.2f} dB<extra></extra>',
                        };
                        layout.xaxis = { title: 'Doppler (m/s)', autorange: true, zeroline: false, showgrid: true };
                        layout.yaxis = { title: 'Range (m)', autorange: true, zeroline: false, showgrid: true };
                    }
                    // --- Handle Time-series or PSD graphs (Scatter plots) ---
                    else if (['raw', 'resp', 'breathPSD', 'heartbeatPSD', 'sensorRaw', 'sensorBreathPSD'].includes(internalGraphTypeKey)) {
                        const rawX = graphData.time || graphData.frequency; // Time for raw/resp/sensorRaw, Frequency for PSD
                        const rawY = graphData.signal || graphData.psd;      // Data for raw/resp/sensorRaw, PSD for PSD

                        const x = flattenArray(rawX);
                        const y = flattenArray(rawY);

                        if (x.length === 0 || y.length === 0 || x.length !== y.length) {
                            console.error(\`createPlotlyTraceForNewWindow: Missing, empty, or mismatched length flattened x/y data for line graph: \${graphData.name} (\${internalGraphTypeKey}).\`);
                            return null;
                        }

                        trace = {
                            ...trace,
                            x: x,
                            y: y,
                            type: 'scatter',
                            mode: 'lines',
                            line: { color: '#0074D9', width: 1.5 },
                            hovertemplate:
                                (internalGraphTypeKey.includes('PSD') ? 'Freq: %{x:.2f} Hz<br>PSD: %{y:.2e}<extra></extra>' : 'Time: %{x:.2f} s<br>Value: %{y:.2f}<extra></extra>')
                        };

                        if (internalGraphTypeKey.includes('PSD')) {
                            layout.xaxis = { title: 'Frequency (Hz)', zeroline: false, showgrid: true };
                            layout.yaxis = { title: 'Power/Magnitude (dB)', type: 'log', zeroline: false, showgrid: true };
                            
                        } else if (['raw', 'resp', 'sensorRaw'].includes(internalGraphTypeKey)) {
                            layout.xaxis = { title: 'Time (s)', zeroline: false, showgrid: true };
                            layout.yaxis = { title: 'Phase Change (rad)' }; // Or 'Amplitude'
                        } else {
                            layout.xaxis = { title: 'X-axis', zeroline: false, showgrid: true };
                            layout.yaxis = { title: 'Y-axis', zeroline: false, showgrid: true };
                        }
                    } else {
                        console.warn(\`createPlotlyTraceForNewWindow: Unsupported graph type or missing data structure for: \${graphData.name} (derived type: \${internalGraphTypeKey})\`);
                        return null;
                    }

                    return { trace, layout };
                }

                                                    // --- Main logic for the new window ---
                                                    const allGraphsString = localStorage.getItem('allAnalysisGraphs');
                                                    const plotsDiv = document.getElementById('plots');

                                                    if (allGraphsString) {
                                                        try {
                                                            const allAnalysisGraphs = JSON.parse(allGraphsString);
                                                            console.log("Retrieved allAnalysisGraphs from localStorage:", allAnalysisGraphs);

                                                            if (allAnalysisGraphs.length === 0) {
                                                                plotsDiv.innerHTML = '<p class="no-data-message">No graph data was found in the analysis results.</p>';
                                                            } else {
                                                                allAnalysisGraphs.forEach((graphData, index) => {
                                                                    const plotlyRenderData = createPlotlyTraceForNewWindow(graphData);
                                                                    if (plotlyRenderData) {
                                                                        const div = document.createElement('div');
                                                                        div.className = "plot-container";
                                                                        div.id = \`graph-\${index}\`; // Unique ID for each plot container
                                                                        plotsDiv.appendChild(div);
                                                                        // Use Plotly.newPlot to render each graph
                                                                        Plotly.newPlot(div, [plotlyRenderData.trace], plotlyRenderData.layout, { responsive: true });
                                                                    } else {
                                                                        // Display a message for graphs that couldn't be rendered
                                                                        const div = document.createElement('div');
                                                                        div.className = "plot-container";
                                                                        div.innerHTML = \`<p class="no-data-message">Could not render: \${graphData.name || 'Unknown Graph'}</p>\`;
                                                                        plotsDiv.appendChild(div);
                                                                    }
                                                                });
                                                            }
                                                        } catch (e) {
                                                            console.error("Error parsing graphs from localStorage:", e);
                                                            plotsDiv.innerHTML = '<p class="no-data-message">Failed to load graph data. Invalid format.</p>';
                                                        } finally {
                                                            // Clean up localStorage after use to avoid stale data on subsequent openings
                                                            localStorage.removeItem('allAnalysisGraphs');
                                                        }
                                                    } else {
                                                        plotsDiv.innerHTML = '<p class="no-data-message">No graph data found in localStorage. Please run an analysis on the main dashboard first.</p>';
                                                    }
                                                </script>
                                            </body>
                                        </html>
                                    `);
                                    graphWindow.document.close(); // Important to close the document stream for proper loading
                                }}
                            >
                                <img src={`${process.env.PUBLIC_URL}/images/logs-light.png`} className="icon" alt="Log" />
                                All Graphs
                            </button>
                        </div>
                    </div>
                    <div
                        id="graphArea"
                        style={{
                        width: '100%',
                        height: '100%',
                        minHeight: '500px',
                        maxHeight: 'calc(100vh - 300px)',
                        overflow: 'hidden',
                        position: 'relative'
                        }}
                        >
                        {/* Conditional rendering for the main graph display */}
                        {mainPlotlyData ? (
                            <Plot
                                key={selectedGraph} // Forces re-render when the selected graph changes
                                data={mainPlotlyData.data}       // Use the dynamically generated trace
                                layout={mainPlotlyData.layout}     // Use the dynamically generated layout
                                useResizeHandler // Enables Plotly to resize with its container
                                style={{ width: "100%", height: "100%" }} // Ensure Plotly fills its container
                            />
                        ) : (
                            <p style={{ textAlign: 'center', marginTop: '50px' }}>
                                {isLoading ? 'Loading graph data...' : 'Select a graph option from the dropdown and click Start.'}
                            </p>
                        )}
                    </div>
                </section>
            </div>
        </div>
    );
};

export default Dashboard;
