import time

def format_time(seconds):
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return hours, minutes, seconds

def log(func):
    def wrapper(*args, **kwargs):
        print(f"\nExecuting {func.__name__}...", flush=True)
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        hours, minutes, seconds = format_time(execution_time)
        print(f"Function '{func.__name__}' completed in {hours:.0f} hours, {minutes:.0f} minutes, {seconds:.4f} seconds.", flush=True)
        return result
    return wrapper
