import React, { useEffect, useState } from 'react';
import { Link } from 'react-router-dom';
import Plot from 'react-plotly.js';
import Plotly from 'plotly.js-dist-min';

function Logs() {
  const [logs, setLogs] = useState([]);
  const [filter, setFilter] = useState('');
  const [date, setDate] = useState('');
  
  // Fetch logs from localStorage on component mount
  useEffect(() => {
    const stored = JSON.parse(localStorage.getItem('logs') || '[]');
    setLogs(stored);
  }, []);
  
  // Filter logs based on search keyword and date
  const filteredLogs = logs.filter(log => {
    const kw = filter.toLowerCase();
    const matchKeyword =
      log.timestamp.toLowerCase().includes(kw) ||
      log.type.toLowerCase().includes(kw) ||
      log.message.toLowerCase().includes(kw);
    const matchDate = date ? log.timestamp.startsWith(date) : true;
    return matchKeyword && matchDate;
  });

  // Function to delete selected logs
  const deleteSelectedLogs = () => {
    const boxes = document.querySelectorAll('.log-checkbox:checked');
    const idxs = Array.from(boxes).map(cb => +cb.dataset.index);
    const newLogs = logs.filter((_, i) => !idxs.includes(i));
    setLogs(newLogs);
    localStorage.setItem('logs', JSON.stringify(newLogs));
  };

  // Function to download log image if it exists
  const downloadImage = dataUrl => {
    const a = document.createElement('a');
    a.href = dataUrl;
    a.download = 'snapshot.png';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
  };

  return (
    <div className="main">
      <div className="content">
        <div className="right-main">
          <section className="log-controls">
            <h2>Logs</h2>
            <div className="log-filter">
              <input
                type="text"
                placeholder="Search logs..."
                value={filter}
                onChange={e => setFilter(e.target.value)}
              />
              <input
                type="date"
                value={date}
                onChange={e => setDate(e.target.value)}
              />
            </div>
            <button
              onClick={deleteSelectedLogs}
              style={{
                marginTop: 10,
                backgroundColor: '#DA291C',
                color: '#fff',
                padding: '10px 16px',
                border: 'none',
                borderRadius: 6,
                cursor: 'pointer'
              }}
            >
              Delete Selected Logs
            </button>
          </section>

          <section className="logs-display">
            <table className="log-table">
              <thead>
                <tr>
                  <th></th>
                  <th>Timestamp</th>
                  <th>Event Type</th>
                  <th>Details</th>
                  <th>Download</th>
                </tr>
              </thead>
              <tbody>
                {filteredLogs.length === 0 ? (
                  <tr>
                    <td colSpan={5}>No logs found.</td>
                  </tr>
                ) : (
                  filteredLogs.map((log, idx) => (
                    <tr key={idx}>
                      <td>
                        <input
                          type="checkbox"
                          className="log-checkbox"
                          data-index={idx}
                        />
                      </td>
                      <td>{log.timestamp}</td>
                      <td className={`type-${log.type.toLowerCase()}`}>
                        {log.type}
                      </td>
                      <td>
                        {log.message}
                        {log.data?.graphType && ` | Graph: ${log.data.graphType}`}
                        {log.data?.uptime && ` | Uptime: ${log.data.uptime}`}
                        {log.data?.scanStatus && ` | Status: ${log.data.scanStatus}`}
                      </td>
                      <td>
                        {log.data?.imageData ? (
                          <button
                            onClick={() => downloadImage(log.data.imageData)}
                            style={{
                              padding: '6px 10px',
                              backgroundColor: '#003087',
                              color: 'white',
                              border: 'none',
                              borderRadius: 4,
                              cursor: 'pointer'
                            }}
                          >
                            Download
                          </button>
                        ) : (
                          <span style={{ color: '#888' }}>—</span>
                        )}
                      </td>
                    </tr>
                  ))
                )}
              </tbody>
            </table>
          </section>
        </div>
      </div>
    </div>
  );
}

export default Logs;
