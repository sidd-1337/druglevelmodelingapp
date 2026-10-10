from __future__ import annotations

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from backend.app.routers import drugs, pk

app = FastAPI(
    title="Drug Level Modeling & Harm Reduction REST API",
    description=(
        "Production backend for pharmacokinetic simulation, SÚKL registry lookup, "
        "and harm reduction dose stacking & safe interval calculations."
    ),
    version="1.0.0",
)

# Enable CORS for Next.js 15 frontend and local tooling
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:3000",
        "http://127.0.0.1:3000",
        "http://localhost:8501",
        "http://127.0.0.1:8501",
        "*",
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Register routers
app.include_router(drugs.router, prefix="/api/v1/drugs", tags=["drugs"])
app.include_router(pk.router, prefix="/api/v1/pk", tags=["pk"])


@app.get("/health", tags=["system"])
def health_check():
    """Health check endpoint for container orchestrators and frontend ping."""
    return {
        "status": "ok",
        "service": "druglevelingapp-backend",
        "version": "1.0.0",
    }


@app.get("/", tags=["system"])
def root():
    """Root info endpoint providing link to interactive API docs."""
    return {
        "message": "Drug Level Modeling & Harm Reduction REST API",
        "docs": "/docs",
        "health": "/health",
    }
