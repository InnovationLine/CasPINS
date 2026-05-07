import os
import json
import logging
from datetime import datetime
from typing import Dict, Any, Optional

class RunLogger:
    """
    Manages structured logging for CasPINS runs.
    Creates detailed, reproducible log files for each analysis module.
    """
    
    def __init__(self, data_dir: str, module_name: str, target_name: str = "unknown"):
        """
        Initialize logger for a specific run.
        
        Args:
            data_dir: Base data directory
            module_name: Name of the module (e.g., 'grna_design', 'primer_design')
            target_name: Name of the target gene/sequence
        """
        self.data_dir = data_dir
        self.module_name = module_name
        self.target_name = target_name.replace(' ', '_').replace('/', '_')
        
        # Create logs directory
        self.logs_dir = os.path.join(data_dir, "logs") if data_dir else os.path.join(os.getcwd(), "logs")
        os.makedirs(self.logs_dir, exist_ok=True)
        
        # Setup log file
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.log_filename = f"{timestamp}_{self.module_name}_{self.target_name}.log"
        self.log_filepath = os.path.join(self.logs_dir, self.log_filename)
        
        # Start run info
        self.run_info = {
            "timestamp_start": datetime.now().isoformat(),
            "module": module_name,
            "target": target_name,
            "parameters": {},
            "status": "running",
            "results_summary": {},
            "errors": []
        }
        
    def log_parameters(self, params: Dict[str, Any]):
        """Record input parameters for the run."""
        self.run_info["parameters"].update(params)
        
    def log_results(self, summary: Dict[str, Any]):
        """Record high-level results/metrics."""
        self.run_info["results_summary"].update(summary)
        
    def log_error(self, error_msg: str):
        """Record an error that occurred during the run."""
        self.run_info["errors"].append({
            "timestamp": datetime.now().isoformat(),
            "message": error_msg
        })
        self.run_info["status"] = "failed"
        
    def finish(self, status: str = "completed") -> str:
        """
        Complete the run log and write to file.
        
        Args:
            status: Final status ('completed', 'failed', etc.)
            
        Returns:
            Path to the written log file
        """
        if self.run_info["status"] != "failed":
            self.run_info["status"] = status
            
        self.run_info["timestamp_end"] = datetime.now().isoformat()
        
        # Calculate duration
        start = datetime.fromisoformat(self.run_info["timestamp_start"])
        end = datetime.fromisoformat(self.run_info["timestamp_end"])
        self.run_info["duration_seconds"] = (end - start).total_seconds()
        
        # Write to file
        with open(self.log_filepath, 'w', encoding='utf-8') as f:
            # Write a human-readable header
            f.write(f"=== CasPINS Run Log ===\n")
            f.write(f"Module: {self.module_name}\n")
            f.write(f"Target: {self.target_name}\n")
            f.write(f"Date: {start.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Status: {self.run_info['status']}\n")
            f.write(f"Duration: {self.run_info['duration_seconds']:.2f}s\n")
            f.write("========================\n\n")
            
            # Write structured JSON payload for reproducibility
            f.write(json.dumps(self.run_info, indent=2))
            
        return self.log_filepath
