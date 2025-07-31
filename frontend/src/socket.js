// src/socket.js
import { io } from 'socket.io-client';

const URL = 'http://localhost:5000';

export const socket = io(URL);

//logging for connection status
socket.on('connect', () => {
  console.log('Connected to Socket.IO server!');
});

socket.on('disconnect', () => {
  console.log('Disconnected from Socket.IO server.');
});

socket.on('connect_error', (err) => {
  console.error('Socket.IO Connection Error:', err.message);
});